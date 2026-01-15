/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MediumSystem.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "MaterialMix.hpp"
#include "MaterialState.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // maximum number of cell densities calculated between two invocations of infoIfElapsed()
    const size_t logProgressChunkSize = 10000;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // This helper class performs spatial sampling of the medium properties within a given cell and
    // for a given medium component. The number of samples to be taken can be specified separately
    // for the density and for all other properties (velocity, magnetic field, metallicity, ...).
    // The main reason for bundling this functionality in its own class/object is to allow caching
    // the random sampling positions with the corresponding sampled density values in a cell.
    // The constructor takes arguments that depend only on the simulation configuration and thus
    // do not vary between cells. This allows consecutively reusing the same object for multiple
    // cells, avoiding memory allocation/deallocation for each cell.
    class PropertySampler
    {
    private:
        vector<Medium*>& _media;
        SpatialGrid* _grid;
        int _numDensitySamples;
        int _numPropertySamples;
        int _numMedia;
        int _cellIndex;
        Position _center;
        vector<Position> _positions;  // indexed on n
        Table<2> _densities;          // indexed on h and n

    public:
        // The constructor allocates the required memory depending on the sampling configuration.
        // The number of samples for each type (density or other properties) can be:
        //   0: the sampling function(s) will never be called
        //   1: "sample" at the cell center only
        //  >1: sample at the specified nr of random points in the cell
        PropertySampler(vector<Medium*>& media, SpatialGrid* grid, int numDensitySamples, int numPropertySamples)
            : _media(media), _grid(grid), _numDensitySamples(numDensitySamples),
              _numPropertySamples(numPropertySamples), _numMedia(_media.size()), _cellIndex(0)
        {
            if (_numDensitySamples > 1 || _numPropertySamples > 1)
            {
                int numSamples = max(_numDensitySamples, _numPropertySamples);
                _positions.resize(numSamples);
                _densities.resize(_numMedia, numSamples);
            }
        }

        // This function should be called with a spatial cell index before calling any of the other
        // functions for that cell. It initializes the required sampling positions and obtains the
        // corresponding sampled density values.
        void prepareForCell(int cellIndex)
        {
            _cellIndex = cellIndex;

            // if we need center "sampling", get the cell center
            if (_numDensitySamples == 1 || _numPropertySamples == 1)
            {
                _center = _grid->centralPositionInCell(_cellIndex);
            }

            // if we need random sampling, draw the positions and get the corresponding densities
            int numSamples = _positions.size();
            for (int n = 0; n != numSamples; ++n)
            {
                _positions[n] = _grid->randomPositionInCell(_cellIndex);
                for (int h = 0; h != _numMedia; ++h) _densities(h, n) = _media[h]->numberDensity(_positions[n]);
            }
        }

        // Returns the central or average density in the cell for the given medium component.
        double density(int h)
        {
            if (_numDensitySamples == 1) return _media[h]->numberDensity(_center);
            double sum = 0.;
            for (int n = 0; n != _numDensitySamples; ++n) sum += _densities(h, n);
            return sum / _numDensitySamples;
        }

        // Returns the central or density-averaged metallicity in the cell for the given medium component.
        double metallicity(int h)
        {
            if (_numPropertySamples == 1) return _media[h]->metallicity(_center);
            double nsum = 0;
            double Zsum = 0.;
            for (int n = 0; n != _numPropertySamples; ++n)
            {
                nsum += _densities(h, n);
                Zsum += _densities(h, n) * _media[h]->metallicity(_positions[n]);
            }
            return nsum > 0. ? Zsum / nsum : 0.;
        }

        // Returns the central or density-averaged temperature in the cell for the given medium component.
        double temperature(int h)
        {
            if (_numPropertySamples == 1) return _media[h]->temperature(_center);
            double nsum = 0;
            double Tsum = 0.;
            for (int n = 0; n != _numPropertySamples; ++n)
            {
                nsum += _densities(h, n);
                Tsum += _densities(h, n) * _media[h]->temperature(_positions[n]);
            }
            return nsum > 0. ? Tsum / nsum : 0.;
        }

        // Returns the central or density-averaged custom parameter values in the cell for the given medium component.
        void parameters(int h, Array& params)
        {
            // always sample at the center to make sure that the output array has the proper size
            _media[h]->parameters(_center, params);
            // if center-sampling is requested, we're done
            if (_numPropertySamples == 1) return;
            // otherwise, clear the output values and accumulate the samples
            params = 0.;
            double nsum = 0;
            Array sample;
            for (int n = 0; n != _numPropertySamples; ++n)
            {
                nsum += _densities(h, n);
                _media[h]->parameters(_positions[n], sample);
                params += _densities(h, n) * sample;
            }
            params /= nsum;
            return;
        }
    };
}

////////////////////////////////////////////////////////////////////

void MediumSystem::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();
    auto log = find<Log>();
    auto parfac = find<ParallelFactory>();
    _config = find<Configuration>();

    _numCells = _grid->numCells();
    if (_numCells < 1) throw FATALERROR("The spatial grid must have at least one cell");
    _numMedia = _media.size();
    size_t allocatedBytes = 0;

    // ----- allocate memory for the medium state -----

    // basic configuration
    _state.initConfiguration(_numCells, _numMedia, 0);

    // common state variables
    vector<StateVariable> variables;
    variables.emplace_back(StateVariable::volume());
    _state.initCommonStateVariables(variables);

    // specific state variables
    for (auto medium : _media) _state.initSpecificStateVariables(medium->mix()->specificStateVariableInfo());

    // finalize
    allocatedBytes += _state.initAllocate() * sizeof(double);

    // ----- allocate memory for the radiation field -----

    if (_config->hasRadiationField())
    {
        _wavelengthGrid = _config->radiationFieldWLG();
        _rf1.resize(_numCells, _wavelengthGrid->numBins());
        allocatedBytes += _rf1.size() * sizeof(double);
    }

    // ----- obtain the material mix pointers -----

    {
        // the material mix pointer is identical for all spatial cells (per medium component)
        _mixv.resize(_numMedia);
        for (int h = 0; h != _numMedia; ++h) _mixv[h] = _media[h]->mix();
    }
    allocatedBytes += _mixv.size() * sizeof(MaterialMix*);

    // cache a list of medium component indices for each material type
    for (int h = 0; h != _numMedia; ++h)
    {
        switch (mix(0, h)->materialType())
        {
            case MaterialMix::MaterialType::Dust: _dust_hv.push_back(h); break;
            case MaterialMix::MaterialType::Gas: _gas_hv.push_back(h); break;
            case MaterialMix::MaterialType::Electrons: _elec_hv.push_back(h); break;
        }
    }

    // ----- inform user about allocated memory -----

    log->info(typeAndName() + " allocated " + StringUtils::toMemSizeString(allocatedBytes) + " of memory");

    // ----- calculate medium properties parallelized on spatial cells -----

    log->info("Calculating medium properties for " + std::to_string(_numCells) + " cells...");
    log->infoSetElapsed(_numCells);
    parfac->parallel()->call(_numCells, [this, log](size_t firstIndex, size_t numIndices) {
        // construct a property sampler to be shared by all cells handled in the loop below
        bool hasProperty = false;
        for (int h = 0; h != _numMedia; ++h)
            if (_media[h]->hasMetallicity() || _media[h]->hasTemperature() || _media[h]->hasParameters())
                hasProperty = true;
        PropertySampler sampler(_media, _grid, _config->numDensitySamples(),
                                hasProperty ? _config->numPropertySamples() : 0);

        // loop over cells
        while (numIndices)
        {
            size_t currentChunkSize = min(logProgressChunkSize, numIndices);
            for (size_t m = firstIndex; m != firstIndex + currentChunkSize; ++m)
            {
                // prepare the sampler for this cell (i.e. store the relevant positions and density samples)
                sampler.prepareForCell(m);

                // volume
                _state.setVolume(m, _grid->volume(m));

                // density: use optional fast-track interface or sample within the cell
                for (int h = 0; h != _numMedia; ++h) _state.setNumberDensity(m, h, sampler.density(h));

                // specific state variables other than density:
                // retrieve value sampled from corresponding medium component
                for (int h = 0; h != _numMedia; ++h)
                {
                    double Z = _media[h]->hasMetallicity() ? sampler.metallicity(h) : -1.;
                    double T = _media[h]->hasTemperature() ? sampler.temperature(h) : -1.;
                    Array params;
                    if (_media[h]->hasParameters()) sampler.parameters(h, params);
                    MaterialState mst(_state, m, h);
                    mix(m, h)->initializeSpecificState(&mst, Z, T, params);
                }
            }
            log->infoIfElapsed("Calculated medium properties: ", currentChunkSize);
            firstIndex += currentChunkSize;
            numIndices -= currentChunkSize;
        }
    });

    // calculate the initial aggregate state, if needed
    _state.calculateAggregate();

    log->info("Done calculating medium properties");
}

////////////////////////////////////////////////////////////////////

int MediumSystem::numMedia() const
{
    return _numMedia;
}

////////////////////////////////////////////////////////////////////

int MediumSystem::numCells() const
{
    return _numCells;
}

////////////////////////////////////////////////////////////////////

const MaterialMix* MediumSystem::mix(int m, int h) const
{
    return _mixPerCell ? _mixv[m * _numMedia + h] : _mixv[h];
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::hasMaterialType(MaterialMix::MaterialType type) const
{
    for (int h = 0; h != _numMedia; ++h)
        if (mix(0, h)->materialType() == type) return true;
    return false;
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::isMaterialType(MaterialMix::MaterialType type, int h) const
{
    return mix(0, h)->materialType() == type;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::dustMassDensity(Position bfr) const
{
    double result = 0.;
    for (int h : _dust_hv) result += _media[h]->massDensity(bfr);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::electronNumberDensity(Position bfr) const
{
    double result = 0.;
    for (int h : _elec_hv) result += _media[h]->numberDensity(bfr);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::gasNumberDensity(Position bfr) const
{
    double result = 0.;
    for (int h : _gas_hv) result += _media[h]->numberDensity(bfr);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::volume(int m) const
{
    return _state.volume(m);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::massDensity(int m) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h) result += massDensity(m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::dustMassDensity(int m) const
{
    double result = 0.;
    for (int h : _dust_hv) result += massDensity(m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::electronNumberDensity(int m) const
{
    double result = 0.;
    for (int h : _elec_hv) result += numberDensity(m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::gasNumberDensity(int m) const
{
    double result = 0.;
    for (int h : _gas_hv) result += numberDensity(m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::numberDensity(int m, int h) const
{
    return _state.numberDensity(m, h);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::massDensity(int m, int h) const
{
    return _state.numberDensity(m, h) * mix(m, h)->mass();
}

////////////////////////////////////////////////////////////////////

double MediumSystem::metallicity(int m, int h) const
{
    return _state.metallicity(m, h);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::temperature(int m, int h) const
{
    return _state.temperature(m, h);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::custom(int m, int h, int i) const
{
    return _state.custom(m, h, i);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityAbs(double lambda, int m, int h) const
{
    MaterialState mst(_state, m, h);
    return mix(m, h)->opacityAbs(lambda, &mst, nullptr);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacitySca(double lambda, int m, int h) const
{
    MaterialState mst(_state, m, h);
    return mix(m, h)->opacitySca(lambda, &mst, nullptr);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m, int h) const
{
    MaterialState mst(_state, m, h);
    return mix(m, h)->opacityExt(lambda, &mst, nullptr);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityAbs(double lambda, int m, MaterialMix::MaterialType type) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h)
        if (mix(0, h)->materialType() == type) result += opacityAbs(lambda, m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m, MaterialMix::MaterialType type) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h)
        if (mix(0, h)->materialType() == type) result += opacityExt(lambda, m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h) result += opacityExt(lambda, m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityAbs(double lambda, int m, const PhotonPacket* pp) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        MaterialState mst(_state, m, h);
        result += mix(m, h)->opacityAbs(lambda, &mst, pp);
    }
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacitySca(double lambda, int m, const PhotonPacket* pp) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        MaterialState mst(_state, m, h);
        result += mix(m, h)->opacitySca(lambda, &mst, pp);
    }
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m, const PhotonPacket* pp) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        MaterialState mst(_state, m, h);
        result += mix(m, h)->opacityExt(lambda, &mst, pp);
    }
    return result;
}

////////////////////////////////////////////////////////////////////
