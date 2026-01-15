/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Configuration.hpp"
#include "FatalError.hpp"
#include "MaterialMix.hpp"
#include "MonteCarloSimulation.hpp"
#include <set>

////////////////////////////////////////////////////////////////////

Configuration::Configuration(SimulationItem* parent)
{
    parent->addChild(this);
}

////////////////////////////////////////////////////////////////////

void Configuration::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // Implementation note: this function is NOT allowed to perform setup on any simulation item in the hierarchy;
    //                      in other words, always use find<XXX>(false) and check for a nullptr result.

    // locate objects that we'll need anyway
    auto sim = find<MonteCarloSimulation>(false);
    if (!sim) throw FATALERROR("Cannot locate a MonteCarloSimulation object in the simulation hierarchy");
    auto ss = find<SourceSystem>(false);
    if (!ss) throw FATALERROR("Cannot locate a SourceSystem object in the simulation hierarchy");

    // ---- set configuration relevant for all simulations, even if there are no media ----

    // retrieve base number of packets
    _numPrimaryPackets = sim->numPackets();

    // retrieve wavelength-related options
    {
        _sourceWavelengthRange.set(ss->minWavelength(), ss->maxWavelength());
        auto is = find<InstrumentSystem>(false);
        if (is) _defaultWavelengthGrid = is->defaultWavelengthGrid();
    }

    // determine the number of media in the simulation hierarchy
    int numMedia = 0;
    auto ms = find<MediumSystem>(false);
    if (ms) numMedia = ms->media().size();  // may be zero
    _hasMedium = (numMedia != 0);

    // if there are no media, we're done
    if (!_hasMedium) return;

    // ---- set configuration relevant for simulations that have media ----

    // retrieve basic photon life-cycle options
    _minWeightReduction = ms->photonPacketOptions()->minWeightReduction();
    _minScattEvents = ms->photonPacketOptions()->minScattEvents();
    _pathLengthBias = ms->photonPacketOptions()->pathLengthBias();

    // retrieve radiation field options
    _hasRadiationField = ms->radiationFieldOptions()->storeRadiationField();
    if (_hasRadiationField)
    {
        _radiationFieldWLG = ms->radiationFieldOptions()->radiationFieldWLG();
    }

    // retrieve media sampling options
    _numDensitySamples = ms->samplingOptions()->numDensitySamples();
    _numPropertySamples = ms->samplingOptions()->numPropertySamples();

    // check for dependencies on extra specific state variables
    bool hasExtraSpecificState = false;
    for (auto medium : ms->media())
        if (medium->mix()->hasExtraSpecificState()) hasExtraSpecificState = true;

    // set the combined medium criteria
    bool _hasConstantSectionMedium = _hasConstantPerceivedWavelength && !hasExtraSpecificState;
    _hasSingleConstantSectionMedium = numMedia == 1 && _hasConstantSectionMedium;
    _hasMultipleConstantSectionMedia = numMedia > 1 && _hasConstantSectionMedium;

    // check for scattering dispersion and for emulating secondary emission
    for (auto medium : ms->media())
    {
        if (medium->mix()->hasScatteringDispersion()) _hasScatteringDispersion = true;
        if (medium->mix()->scatteringEmulatesSecondaryEmission()) _scatteringEmulatesSecondaryEmission = true;
    }
    _needIndividualPeelOff = _hasScatteringDispersion | _scatteringEmulatesSecondaryEmission;
}

////////////////////////////////////////////////////////////////////

void Configuration::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    // Implementation note: this function should perform pure logging and is NOT allowed to perform setup
    //                      on any simulation item in the hierarchy except for the logger.
    auto log = find<Log>();

    // --- log wavelength regime, simulation mode, and media characteristics  ---

    string regime = _oligochromatic ? "Oligo" : "Pan";
    log->info("  " + regime + "chromatic wavelength regime");
    string medium = _hasMedium ? "With" : "No";
    log->info("  " + medium + " transfer medium");
}

////////////////////////////////////////////////////////////////////

namespace
{
    // This function extends the specified wavelength range with the range of the specified wavelength grid
    void extendForWavelengthGrid(Range& range, WavelengthGrid* wavelengthGrid)
    {
        // we explicitly call setup() on wavelength grids before accessing them
        // because this function may be called early during simulation setup
        wavelengthGrid->setup();
        range.extend(wavelengthGrid->wavelengthRange());
    }
}

////////////////////////////////////////////////////////////////////

Range Configuration::simulationWavelengthRange() const
{
    // include primary and secondary source ranges
    Range range = _sourceWavelengthRange;

    // include radiation field wavelength grid (because dust properties are pre-calculated on these wavelengths)
    if (_hasRadiationField)
    {
        extendForWavelengthGrid(range, _radiationFieldWLG);
        // the calculation of the Planck-integrated absorption cross sections needs this wavelength range;
        // see the precalculate() function in EquilibriumDustEmissionCalculator and StochasticDustEmissionCalculator
        range.extend(Range(0.09e-6, 2000e-6));
    }

    // include default instrument wavelength grid
    if (_defaultWavelengthGrid) extendForWavelengthGrid(range, _defaultWavelengthGrid);

    // include instrument-specific wavelength grids
    auto is = find<InstrumentSystem>(false);
    for (auto ins : is->instruments())
    {
        if (ins->wavelengthGrid()) extendForWavelengthGrid(range, ins->wavelengthGrid());
    }

    // extend the final range with a narrow margin for round-offs
    range.extendWithRedshift(1. / 100.);
    return range;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // This function adds the characteristic wavelengths of the specified grid to the specified set of wavelengths
    void addForWavelengthGrid(std::set<double>& wavelengths, WavelengthGrid* wavelengthGrid)
    {
        // we explicitly call setup() on wavelength grids before accessing them
        // because this function may be called early during simulation setup
        wavelengthGrid->setup();
        int n = wavelengthGrid->numBins();
        for (int ell = 0; ell != n; ++ell) wavelengths.insert(wavelengthGrid->wavelength(ell));
    }
}

////////////////////////////////////////////////////////////////////

vector<double> Configuration::simulationWavelengths() const
{
    std::set<double> wavelengths;

    // include radiation field and dust emission wavelength grids
    if (_hasRadiationField) addForWavelengthGrid(wavelengths, _radiationFieldWLG);

    // include default instrument wavelength grid
    if (_defaultWavelengthGrid) addForWavelengthGrid(wavelengths, _defaultWavelengthGrid);

    // include instrument-specific wavelength grids
    auto is = find<InstrumentSystem>(false);
    for (auto ins : is->instruments())
    {
        if (ins->wavelengthGrid()) addForWavelengthGrid(wavelengths, ins->wavelengthGrid());
    }

    return vector<double>(wavelengths.begin(), wavelengths.end());
}

////////////////////////////////////////////////////////////////////

WavelengthGrid* Configuration::wavelengthGrid(WavelengthGrid* localWavelengthGrid) const
{
    auto result = localWavelengthGrid && !_oligochromatic ? localWavelengthGrid : _defaultWavelengthGrid;
    if (!result) throw FATALERROR("Cannot find a wavelength grid for instrument or probe");
    return result;
}

////////////////////////////////////////////////////////////////////
