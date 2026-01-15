/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SimulationItemRegistry.hpp"
#include "ItemRegistry.hpp"
#include "SkirtUnitDef.hpp"

// ---> add new items below in alphabetical order

#include "CartesianSpatialGrid.hpp"
#include "CubicSplineSmoothingKernel.hpp"
#include "ExtragalacticUnits.hpp"
#include "FileSED.hpp"
#include "InstrumentSystem.hpp"
#include "LinMesh.hpp"
#include "LogWavelengthGrid.hpp"
#include "MeanInterstellarDustMix.hpp"
#include "MediumSystem.hpp"
#include "MonteCarloSimulation.hpp"
#include "ParticleMedium.hpp"
#include "ParticleSource.hpp"
#include "ProbeSystem.hpp"
#include "Random.hpp"
#include "SEDInstrument.hpp"
#include "SIUnits.hpp"
#include "SourceSystem.hpp"
#include "SpatialGrid.hpp"
#include "StellarUnits.hpp"

////////////////////////////////////////////////////////////////////

SimulationItemRegistry::SimulationItemRegistry(string version, string format)
{
    // start a new schema
    ItemRegistry::beginSchema("SKIRT", "a SKIRT parameter file", version, "ski", "skirt-simulation-hierarchy",
                              "MonteCarloSimulation", format, "http://www.skirt.ugent.be/skirt9");

    // add the SKIRT unit definitions
    ItemRegistry::addUnitDef<SkirtUnitDef>();

    // add the SKIRT simulation items
    ItemRegistry::add<SimulationItem>();

    // ---> add new items in the order you want them to appear in choice lists for the user

    // basic building blocks
    ItemRegistry::add<Simulation>();
    ItemRegistry::add<Random>();
    ItemRegistry::add<Units>();
    ItemRegistry::add<SIUnits>();
    ItemRegistry::add<StellarUnits>();
    ItemRegistry::add<ExtragalacticUnits>();

    // source system and sources
    ItemRegistry::add<SourceSystem>();
    ItemRegistry::add<Source>();
    ItemRegistry::add<ImportedSource>();
    ItemRegistry::add<ParticleSource>();

    // luminosity normalizations

    // SEDs
    ItemRegistry::add<SED>();
    ItemRegistry::add<ContSED>();
    ItemRegistry::add<FileSED>();

    // SED families
    ItemRegistry::add<SEDFamily>();
    // GPU-SKIRT: choose one here
    // ItemRegistry::add<BlackBodySEDFamily>();
    // ItemRegistry::add<CastelliKuruczSEDFamily>();
    // ItemRegistry::add<BruzualCharlotSEDFamily>();
    // ItemRegistry::add<MarastonSEDFamily>();
    // ItemRegistry::add<Starburst99SEDFamily>();
    // ItemRegistry::add<FSPSSEDFamily>();
    // ItemRegistry::add<BpassSEDFamily>();
    // ItemRegistry::add<FileSSPSEDFamily>();
    // ItemRegistry::add<FileIndexedSEDFamily>();
    // ItemRegistry::add<MappingsSEDFamily>();
    // ItemRegistry::add<ToddlersSEDFamily>();
    // ItemRegistry::add<SpinFlipSEDFamily>();
    // ItemRegistry::add<LyaGaussianSEDFamily>();
    // ItemRegistry::add<LyaDoublePeakedSEDFamily>();
    // ItemRegistry::add<LyaSEDFamilyDecorator>();

    // wavelength distributions
    ItemRegistry::add<WavelengthDistribution>();
    // GPU-SKIRT: choose one here
    // ItemRegistry::add<DefaultWavelengthDistribution>();
    // ItemRegistry::add<RangeWavelengthDistribution>();
    // ItemRegistry::add<LinWavelengthDistribution>();
    // ItemRegistry::add<LogWavelengthDistribution>();
    // ItemRegistry::add<TabulatedWavelengthDistribution>();
    // ItemRegistry::add<FileWavelengthDistribution>();
    // ItemRegistry::add<ListWavelengthDistribution>();
    // ItemRegistry::add<DiscreteWavelengthDistribution>();

    // bands

    // angular distributions

    // polarization profiles

    // geometries

    // geometry decorators

    // smoothing kernels
    ItemRegistry::add<SmoothingKernel>();
    ItemRegistry::add<CubicSplineSmoothingKernel>();

    // vector fields

    // spatial grids
    ItemRegistry::add<SpatialGrid>();
    ItemRegistry::add<BoxSpatialGrid>();
    ItemRegistry::add<CartesianSpatialGrid>();

    // spatial grid policies

    // one-dimensional meshes for spatial grids
    ItemRegistry::add<Mesh>();
    ItemRegistry::add<LinMesh>();

    // medium system and media
    ItemRegistry::add<MediumSystem>();
    ItemRegistry::add<Medium>();
    ItemRegistry::add<ImportedMedium>();
    ItemRegistry::add<ParticleMedium>();

    // medium system options
    ItemRegistry::add<PhotonPacketOptions>();
    ItemRegistry::add<RadiationFieldOptions>();
    ItemRegistry::add<SamplingOptions>();

    // material normalizations
    // GPU-SKIRT: choose one here
    // ItemRegistry::add<MassMaterialNormalization>();
    // ItemRegistry::add<NumberMaterialNormalization>();
    // ItemRegistry::add<AxisMaterialNormalization>();
    // ItemRegistry::add<OpticalDepthMaterialNormalization>();
    // ItemRegistry::add<MassColumnMaterialNormalization>();
    // ItemRegistry::add<NumberColumnMaterialNormalization>();

    // material mixes
    ItemRegistry::add<MaterialMix>();
    ItemRegistry::add<DustMix>();

    ItemRegistry::add<SingleGrainDustMix>();
    ItemRegistry::add<MeanInterstellarDustMix>();

    // material mix families

    // grain population
    // GPU-SKIRT: choose one here
    // ItemRegistry::add<GrainPopulation>();

    // grain size distributions
    // GPU-SKIRT: choose one here
    // ItemRegistry::add<GrainSizeDistribution>();
    // ItemRegistry::add<RangeGrainSizeDistribution>();
    // ItemRegistry::add<PowerLawGrainSizeDistribution>();
    // ItemRegistry::add<ModifiedPowerLawGrainSizeDistribution>();
    // ItemRegistry::add<LogNormalGrainSizeDistribution>();
    // ItemRegistry::add<ModifiedLogNormalGrainSizeDistribution>();
    // ItemRegistry::add<SingleGrainSizeDistribution>();
    // ItemRegistry::add<ZubkoSilicateGrainSizeDistribution>();
    // ItemRegistry::add<ZubkoGraphiteGrainSizeDistribution>();
    // ItemRegistry::add<ZubkoPAHGrainSizeDistribution>();
    // ItemRegistry::add<HirashitaLogNormalGrainSizeDistribution>();
    // ItemRegistry::add<FileGrainSizeDistribution>();
    // ItemRegistry::add<ListGrainSizeDistribution>();

    // grain compositions
    // GPU-SKIRT: choose one here
    // ItemRegistry::add<GrainComposition>();
    // ItemRegistry::add<DraineSilicateGrainComposition>();
    // ItemRegistry::add<DraineGraphiteGrainComposition>();
    // ItemRegistry::add<DraineNeutralPAHGrainComposition>();
    // ItemRegistry::add<DraineIonizedPAHGrainComposition>();
    // ItemRegistry::add<MieSilicateGrainComposition>();
    // ItemRegistry::add<MinSilicateGrainComposition>();
    // ItemRegistry::add<PolarizedSilicateGrainComposition>();
    // ItemRegistry::add<PolarizedGraphiteGrainComposition>();
    // ItemRegistry::add<CrystalEnstatiteGrainComposition>();
    // ItemRegistry::add<CrystalForsteriteGrainComposition>();
    // ItemRegistry::add<DustEmGrainComposition>();
    // ItemRegistry::add<BegemannPorousAluminaGrainComposition>();
    // ItemRegistry::add<HofmeisterPericlaseGrainComposition>();
    // ItemRegistry::add<DorschnerOlivineGrainComposition>();
    // ItemRegistry::add<TrustSilicateGrainComposition>();
    // ItemRegistry::add<TrustGraphiteGrainComposition>();
    // ItemRegistry::add<TrustNeutralPAHGrainComposition>();
    // ItemRegistry::add<SpheroidalSilicateGrainComposition>();
    // ItemRegistry::add<SpheroidalGraphiteGrainComposition>();

    // spatial cell libraries

    // dynamic medium state recipes

    // wavelength grids
    ItemRegistry::add<WavelengthGrid>();
    ItemRegistry::add<DisjointWavelengthGrid>();
    ItemRegistry::add<LogWavelengthGrid>();

    // instrument system and instruments
    ItemRegistry::add<InstrumentSystem>();
    ItemRegistry::add<Instrument>();
    ItemRegistry::add<DistantInstrument>();
    ItemRegistry::add<SEDInstrument>();

    // all-sky projections

    // probe system and probes
    ItemRegistry::add<ProbeSystem>();
    ItemRegistry::add<Probe>();
    //   .. convergence
    //   .. source
    //   .. spatial grid
    //   .. properties
    //   .. wavelength grid
    //   .. specialty
    //   .. imported source
    //   .. imported medium

    // forms

    // Monte Carlo simulations
    ItemRegistry::add<MonteCarloSimulation>();
}

////////////////////////////////////////////////////////////////////

const SchemaDef* SimulationItemRegistry::getSchemaDef()
{
    return ItemRegistry::getSchemaDef("SKIRT");
}

////////////////////////////////////////////////////////////////////

SimulationItemRegistry::~SimulationItemRegistry()
{
    ItemRegistry::finalize();
}

////////////////////////////////////////////////////////////////////
