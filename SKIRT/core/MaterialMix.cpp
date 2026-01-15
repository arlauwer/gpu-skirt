/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MaterialMix.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void MaterialMix::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    _random = find<Random>();
    _config = find<Configuration>();
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasExtraSpecificState() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasScatteringDispersion() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::scatteringEmulatesSecondaryEmission() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> MaterialMix::parameterInfo() const
{
    return vector<SnapshotParameter>();
}

////////////////////////////////////////////////////////////////////

void MaterialMix::initializeSpecificState(MaterialState* /*state*/, double /*metallicity*/, double /*temperature*/,
                                          const Array& /*params*/) const
{}

////////////////////////////////////////////////////////////////////

UpdateStatus MaterialMix::updateSpecificState(MaterialState* /*state*/, const Array& /*Jv*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::isSpecificStateConverged(int /*numCells*/, int /*numUpdated*/, int /*numNotConverged*/,
                                           MaterialState* /*currentAggregate*/,
                                           MaterialState* /*previousAggregate*/) const
{
    return true;
}

////////////////////////////////////////////////////////////////////

double MaterialMix::asymmpar(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::peeloffScattering(double& /*I*/, double& /*Q*/, double& /*U*/, double& /*V*/, double& /*lambda*/,
                                    Direction /*bfkobs*/, Direction /*bfky*/, const MaterialState* /*state*/,
                                    const PhotonPacket* /*pp*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

void MaterialMix::performScattering(double /*lambda*/, const MaterialState* /*state*/, PhotonPacket* /*pp*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////
