/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MediumState.hpp"
#include "FatalError.hpp"

//////////////////////////////////////////////////////////////////////

void MediumState::initConfiguration(int numCells, int numMedia, int numAggregateCells)
{
    _numCells = numCells;
    _numMedia = numMedia;
    _numAggregateCells = numAggregateCells;

    _off_dens.resize(_numMedia);
    _off_meta.resize(_numMedia);
    _off_temp.resize(_numMedia);
    _off_cust.resize(_numMedia);
}

//////////////////////////////////////////////////////////////////////

void MediumState::initCommonStateVariables(const vector<StateVariable>& variables)
{
    for (const StateVariable& variable : variables)
    {
        switch (variable.identifier())
        {
            case StateVariable::Identifier::Volume: _off_volu = _nextOffset++; break;
            case StateVariable::Identifier::BulkVelocity:
                _off_velo = _nextOffset;
                _nextOffset += 3;
                break;
            case StateVariable::Identifier::MagneticField:
                _off_mfld = _nextOffset;
                _nextOffset += 3;
                break;
            case StateVariable::Identifier::NumberDensity:
            case StateVariable::Identifier::Metallicity:
            case StateVariable::Identifier::Temperature:
            case StateVariable::Identifier::Custom:
                throw FATALERROR("Requesting common state variable of unsupported type");
        }
    }
}

//////////////////////////////////////////////////////////////////////

void MediumState::initSpecificStateVariables(const vector<StateVariable>& variables)
{
    for (const StateVariable& variable : variables)
    {
        switch (variable.identifier())
        {
            case StateVariable::Identifier::NumberDensity:
                if (_numAggregateCells) _densityOffsets.push_back(_nextOffset);
                _off_dens[_nextComponent] = _nextOffset++;
                break;
            case StateVariable::Identifier::Metallicity: _off_meta[_nextComponent] = _nextOffset++; break;
            case StateVariable::Identifier::Temperature: _off_temp[_nextComponent] = _nextOffset++; break;
            case StateVariable::Identifier::Custom:
                if (variable.customIndex() == 0) _off_cust[_nextComponent] = _nextOffset;
                if (_numAggregateCells && variable.quantity() == "numbervolumedensity")
                    _densityOffsets.push_back(_nextOffset);
                _nextOffset++;
                break;
            case StateVariable::Identifier::Volume:
            case StateVariable::Identifier::BulkVelocity:
            case StateVariable::Identifier::MagneticField:
                throw FATALERROR("Requesting specific state variable of unsupported type");
        }
    }
    _nextComponent++;
}

//////////////////////////////////////////////////////////////////////

size_t MediumState::initAllocate()
{
    if (_nextComponent != _numMedia) throw FATALERROR("Failed to request state variables for all medium components");
    _numVars = _nextOffset;

    size_t numAlloc = static_cast<size_t>(_numVars) * static_cast<size_t>(_numCells + _numAggregateCells);
    _data.resize(numAlloc);
    return numAlloc;
}

//////////////////////////////////////////////////////////////////////

std::pair<int, int> MediumState::synchronize(const vector<UpdateStatus>& cellFlags)
{
    int numUpdated = 0;
    int numNotConverged = 0;

    for (int m = 0; m != _numCells; ++m)
    {
        if (cellFlags[m].isUpdated()) numUpdated++;
        if (!cellFlags[m].isConverged()) numNotConverged++;
    }

    return std::make_pair(numUpdated, numNotConverged);
}

//////////////////////////////////////////////////////////////////////

void MediumState::calculateAggregate()
{
    // if aggregation is requested
    if (_numAggregateCells)
    {
        // clear the variables in the current aggregate state
        for (int i = 0; i != _numVars; ++i) _data[_numVars * _numCells + i] = 0.;

        // calculate the current aggregate state
        for (int m = 0; m != _numCells; ++m)
        {
            // cell volume
            double volume = _data[_numVars * m + _off_volu];
            _data[_numVars * _numCells + _off_volu] += volume;

            // variables of type number volume density
            for (int i : _densityOffsets) _data[_numVars * _numCells + i] += _data[_numVars * m + i] * volume;
        }
    }
}

//////////////////////////////////////////////////////////////////////

void MediumState::pushAggregate()
{
    // if aggregation is requested
    if (_numAggregateCells)
    {
        // shift the previous aggregate states to make room
        for (int m = _numCells + _numAggregateCells - 1; m != _numCells; --m)
        {
            for (int i = 0; i != _numVars; ++i) _data[_numVars * m + i] = _data[_numVars * (m - 1) + i];
        }
    }
}

//////////////////////////////////////////////////////////////////////

void MediumState::setVolume(int m, double value)
{
    _data[_numVars * m + _off_volu] = value;
}

//////////////////////////////////////////////////////////////////////

void MediumState::setBulkVelocity(int m, Vec value)
{
    int i = _numVars * m + _off_velo;
    _data[i] = value.x();
    _data[i + 1] = value.y();
    _data[i + 2] = value.z();
}

//////////////////////////////////////////////////////////////////////

void MediumState::setMagneticField(int m, Vec value)
{
    int i = _numVars * m + _off_mfld;
    _data[i] = value.x();
    _data[i + 1] = value.y();
    _data[i + 2] = value.z();
}

//////////////////////////////////////////////////////////////////////

void MediumState::setNumberDensity(int m, int h, double value)
{
    _data[_numVars * m + _off_dens[h]] = value;
}

//////////////////////////////////////////////////////////////////////

void MediumState::setMetallicity(int m, int h, double value)
{
    _data[_numVars * m + _off_meta[h]] = value;
}

//////////////////////////////////////////////////////////////////////

void MediumState::setTemperature(int m, int h, double value)
{
    _data[_numVars * m + _off_temp[h]] = value;
}

//////////////////////////////////////////////////////////////////////

void MediumState::setCustom(int m, int h, int i, double value)
{
    _data[_numVars * m + _off_cust[h] + i] = value;
}

//////////////////////////////////////////////////////////////////////
