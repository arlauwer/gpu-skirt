/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustMix.hpp"
#include "Configuration.hpp"
#include "Log.hpp"
#include "MaterialState.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

// built-in constants determining the resolution for discretizing the scattering angles
namespace
{
    // theta from 0 to pi, index t
    constexpr int numTheta = 361;
    constexpr int maxTheta = numTheta - 1;
    constexpr double deltaTheta = M_PI / maxTheta;

    // phi from 0 to 2 pi, index f
    constexpr int numPhi = 361;
    constexpr int maxPhi = numPhi - 1;
    constexpr double deltaPhi = 2 * M_PI / (maxPhi - 1);
}

// hard-coded parameters configuring extreme forward/backward Henyey-Greenstein scattering
namespace
{
    // calculating the value of and sampling from the HG phase function is numerically unstable for abs(g) > gmax
    constexpr double gmax = 0.999999;
}

////////////////////////////////////////////////////////////////////

void DustMix::setupSelfAfter()
{
    MaterialMix::setupSelfAfter();
    auto config = find<Configuration>();

    // construct a wavelength grid for sampling dust properties containing two (merged) sets of grid points:
    //  - a fine grid (in log space) covering the entire wavelength range of the simulation
    //  - all specific wavelengths mentioned in the configuration of the simulation (grids, normalizations, ...)
    //    ensuring that the dust properties are interpolated at exactly these wavelengths (e.g. for benchmarks)
    // we first gather all the wavelength points, in arbitrary order, and then sort them
    vector<double> wavelengths;
    wavelengths.reserve(10000);

    // add a fine grid covering the wavelength range of the simulation in log space;
    // use integer multiples as logarithmic grid points so that the grid is stable for changing wavelength ranges
    const double numWavelengthsPerDex = 1000;  // declared as double to avoid casting, but should be an integer
    Range range = config->simulationWavelengthRange();
    int minLambdaSerial = std::floor(numWavelengthsPerDex * log10(range.min()));
    int maxLambdaSerial = std::ceil(numWavelengthsPerDex * log10(range.max()));
    for (int k = minLambdaSerial; k <= maxLambdaSerial; ++k) wavelengths.push_back(pow(10., k / numWavelengthsPerDex));

    // add the wavelengths mentioned in the configuration of the simulation
    for (double lambda : config->simulationWavelengths()) wavelengths.push_back(lambda);

    // sort the wavelengths and remove duplicates
    NR::unique(wavelengths);

    // remove wavelengths beyond 10 cm, if any, because those dust cross sections should be assumed to be zero
    const double dm = 0.1;
    bool radioCutoff = wavelengths.back() > dm;
    if (radioCutoff)
    {
        wavelengths.resize(NR::locate(wavelengths, dm) + 1);
        if (wavelengths.empty() || wavelengths.back() != dm) wavelengths.push_back(dm);
        wavelengths.push_back(dm * 1.001);  // add a dummy border value which is never used
    }

    // remember the wavelength range
    _required.set(wavelengths.front(), wavelengths.back());

    // copy the wavelengths into a temporary (local) array that will be used to sample the dust properties
    Array lambdav = NR::array(wavelengths);
    int numLambda = lambdav.size();

    // derive a wavelength grid that will be used for converting a wavelength to an index in the above array;
    // the grid points are shifted to the left of the actual sample points to approximate rounding
    _lambdav.resize(numLambda);
    _lambdav[0] = lambdav[0];
    for (int ell = 1; ell != numLambda; ++ell)
    {
        _lambdav[ell] = sqrt(lambdav[ell] * lambdav[ell - 1]);
    }

    // get the scattering mode advertised by this dust mix
    auto mode = scatteringMode();

    // if needed, build a scattering angle grid
    if (mode == ScatteringMode::MaterialPhaseFunction || mode == ScatteringMode::SphericalPolarization
        || mode == ScatteringMode::SpheroidalPolarization)
    {
        _thetav.resize(numTheta);
        for (int t = 0; t != numTheta; ++t) _thetav[t] = t * deltaTheta;
    }

    // resize the optical property arrays and tables as needed
    _sigmaabsv.resize(numLambda);
    _sigmascav.resize(numLambda);
    _sigmaextv.resize(numLambda);
    _asymmparv.resize(numLambda);
    if (mode == ScatteringMode::MaterialPhaseFunction || mode == ScatteringMode::SphericalPolarization
        || mode == ScatteringMode::SpheroidalPolarization)
    {
        _S11vv.resize(numLambda, numTheta);
        if (mode == ScatteringMode::SphericalPolarization || mode == ScatteringMode::SpheroidalPolarization)
        {
            _S12vv.resize(numLambda, numTheta);
            _S33vv.resize(numLambda, numTheta);
            _S34vv.resize(numLambda, numTheta);
        }
        if (mode == ScatteringMode::SpheroidalPolarization)
        {
            _sigmaabsvv.resize(numLambda, numTheta);
            _sigmaabspolvv.resize(numLambda, numTheta);
        }
    }

    // obtain the optical properties from the subclass
    _mu = getOpticalProperties(lambdav, _thetav, _sigmaabsv, _sigmascav, _asymmparv, _S11vv, _S12vv, _S33vv, _S34vv,
                               _sigmaabsvv, _sigmaabspolvv);

    // ensure that values for the asymmetry parameter g are not too close to 1 or -1
    int lastAdjustedEll = -1;  // index of longest wavelength for which adjustment has been made
    for (int ell = 0; ell != numLambda; ++ell)
    {
        if (abs(_asymmparv[ell]) > gmax)
        {
            _asymmparv[ell] = std::copysign(gmax, _asymmparv[ell]);
            lastAdjustedEll = ell;
        }
    }
    if (lastAdjustedEll >= 0)
    {
        auto units = find<Units>();
        auto log = find<Log>();
        log->warning(type() + " extreme forward/backward scattering asymmetry parameter values reduced to "
                     + StringUtils::toString(gmax));
        log->warning("  For some or all wavelengths short of "
                     + StringUtils::toString(units->owavelength(lambdav[lastAdjustedEll]), 'g', 3) + " "
                     + units->uwavelength());
    }

    // calculate derived basic optical properties
    for (int ell = 0; ell != numLambda; ++ell)
    {
        _sigmaextv[ell] = _sigmaabsv[ell] + _sigmascav[ell];
    }

    // precalculate discretizations related to the scattering angles as needed
    if (mode == ScatteringMode::MaterialPhaseFunction || mode == ScatteringMode::SphericalPolarization
        || mode == ScatteringMode::SpheroidalPolarization)
    {
        // create a table with the normalized cumulative distribution of theta for each wavelength
        _thetaXvv.resize(numLambda, 0);
        for (int ell = 0; ell != numLambda; ++ell)
        {
            NR::cdf(_thetaXvv[ell], maxTheta, [this, ell](int t) { return _S11vv(ell, t + 1) * sin(_thetav[t + 1]); });
        }

        // create a table with the phase function normalization factor for each wavelength
        _pfnormv.resize(numLambda);
        for (int ell = 0; ell != numLambda; ++ell)
        {
            double sum = 0.;
            for (int t = 0; t != numTheta; ++t)
            {
                sum += _S11vv(ell, t) * sin(_thetav[t]) * deltaTheta;
            }
            _pfnormv[ell] = 2.0 / sum;
        }

        // create tables listing phi, phi/(2 pi), sin(2 phi) and 1-cos(2 phi) for each phi index
        if (mode == ScatteringMode::SphericalPolarization || mode == ScatteringMode::SpheroidalPolarization)
        {
            _phiv.resize(numPhi);
            _phi1v.resize(numPhi);
            _phisv.resize(numPhi);
            _phicv.resize(numPhi);
            for (int f = 0; f != numPhi; ++f)
            {
                double phi = f * deltaPhi;
                _phiv[f] = phi;
                _phi1v[f] = phi / (2 * M_PI);
                _phisv[f] = sin(2 * phi);
                _phicv[f] = 1 - cos(2 * phi);
            }
        }
    }

    // if the wavelength range was cut off, suppress all cross sections beyond the cutoff wavelength
    // by setting the last two values to zero (the penultimate item is for the cutoff wavelength)
    if (radioCutoff)
    {
        _sigmaabsv[numLambda - 2] = _sigmaabsv[numLambda - 1] = 0.;
        _sigmascav[numLambda - 2] = _sigmascav[numLambda - 1] = 0.;
        _sigmaextv[numLambda - 2] = _sigmaextv[numLambda - 1] = 0.;
    }

    // give the subclass a chance to obtain additional precalculated information
    size_t allocatedBytes = initializeExtraProperties(lambdav);

    // calculate and log allocated memory size
    size_t allocatedSize = 0;
    allocatedSize += _thetav.size();
    allocatedSize += _sigmaabsv.size();
    allocatedSize += _sigmascav.size();
    allocatedSize += _sigmaextv.size();
    allocatedSize += _asymmparv.size();
    allocatedSize += _S11vv.size();
    allocatedSize += _S12vv.size();
    allocatedSize += _S33vv.size();
    allocatedSize += _S34vv.size();
    allocatedSize += _thetaXvv.size();
    allocatedSize += _pfnormv.size();
    allocatedSize += _phiv.size();
    allocatedSize += _phi1v.size();
    allocatedSize += _phisv.size();
    allocatedSize += _phicv.size();
    allocatedSize += _sigmaabsvv.size();
    allocatedSize += _sigmaabspolvv.size();

    allocatedBytes += allocatedSize * sizeof(double);
    find<Log>()->info(type() + " allocated " + StringUtils::toMemSizeString(allocatedBytes) + " of memory");
}

////////////////////////////////////////////////////////////////////

size_t DustMix::initializeExtraProperties(const Array& /*lambdav*/)
{
    return 0;
}

////////////////////////////////////////////////////////////////////

void DustMix::informAvailableWavelengthRange(Range available)
{
    const double fuzzy = 0.011;  // slightly more than 1% because Configuration class extends the range by 1%
    if (!available.containsFuzzy(_required.min(), fuzzy) || !available.containsFuzzy(_required.max(), fuzzy))
    {
        auto units = find<Units>();
        auto outstring = [units](double wavelength) {
            return StringUtils::toString(units->owavelength(wavelength), 'g', 3) + " " + units->uwavelength();
        };

        auto log = find<Log>();
        log->warning(type() + " dust properties are not available for full simulation wavelength range");
        log->warning("  Required: " + outstring(_required.min()) + " - " + outstring(_required.max()));
        log->warning("  Available: " + outstring(available.min()) + " - " + outstring(available.max()));
    }
}

////////////////////////////////////////////////////////////////////

int DustMix::indexForLambda(double lambda) const
{
    return NR::locateClip(_lambdav, lambda);
}

////////////////////////////////////////////////////////////////////

int DustMix::indexForTheta(double theta) const
{
    int t = 0.5 + theta / deltaTheta;
    return max(0, min(t, maxTheta));
}

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType DustMix::materialType() const
{
    return MaterialType::Dust;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> DustMix::specificStateVariableInfo() const
{
    return vector<StateVariable>{StateVariable::numberDensity()};
}

////////////////////////////////////////////////////////////////////

double DustMix::mass() const
{
    return _mu;
}

////////////////////////////////////////////////////////////////////

double DustMix::sectionAbs(double lambda) const
{
    return _sigmaabsv[indexForLambda(lambda)];
}

////////////////////////////////////////////////////////////////////

double DustMix::sectionSca(double lambda) const
{
    return _sigmascav[indexForLambda(lambda)];
}

////////////////////////////////////////////////////////////////////

double DustMix::sectionExt(double lambda) const
{
    return _sigmaextv[indexForLambda(lambda)];
}

////////////////////////////////////////////////////////////////////

double DustMix::asymmpar(double lambda) const
{
    return _asymmparv[indexForLambda(lambda)];
}

////////////////////////////////////////////////////////////////////

double DustMix::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double n = state->numberDensity();
    return n > 0. ? n * _sigmaabsv[indexForLambda(lambda)] : 0.;
}

////////////////////////////////////////////////////////////////////

double DustMix::opacitySca(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double n = state->numberDensity();
    return n > 0. ? n * _sigmascav[indexForLambda(lambda)] : 0.;
}

////////////////////////////////////////////////////////////////////

double DustMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double n = state->numberDensity();
    return n > 0. ? n * _sigmaextv[indexForLambda(lambda)] : 0.;
}

////////////////////////////////////////////////////////////////////

DustMix::ScatteringMode DustMix::scatteringMode() const
{
    return DustMix::ScatteringMode::HenyeyGreenstein;
}

////////////////////////////////////////////////////////////////////

double DustMix::phaseFunctionValueForCosine(double lambda, double costheta) const
{
    int ell = indexForLambda(lambda);
    int t = indexForTheta(acos(costheta));
    return _pfnormv[ell] * _S11vv(ell, t);
}

////////////////////////////////////////////////////////////////////

double DustMix::generateCosineFromPhaseFunction(double lambda) const
{
    return cos(random()->cdfLinLin(_thetav, _thetaXvv[indexForLambda(lambda)]));
}

////////////////////////////////////////////////////////////////////
