/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SEDInstrument.hpp"

////////////////////////////////////////////////////////////////////

void SEDInstrument::setupSelfBefore()
{
    DistantInstrument::setupSelfBefore();

    // configure flux recorder
    // GPU-SKIRT

    // precalculate information needed by detect() function
    _radius2 = radius() * radius();
    _costheta = cos(inclination());
    _sintheta = sin(inclination());
    _cosphi = cos(azimuth());
    _sinphi = sin(azimuth());
}

////////////////////////////////////////////////////////////////////
