/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MONTECARLOSIMULATION_HPP
#define MONTECARLOSIMULATION_HPP

#include "Configuration.hpp"
#include "InstrumentSystem.hpp"
#include "MediumSystem.hpp"
#include "ProbeSystem.hpp"
#include "Simulation.hpp"
#include "SourceSystem.hpp"
class SecondarySourceSystem;

//////////////////////////////////////////////////////////////////////

/** The MonteCarloSimulation class is the top-level class describing a SKIRT simulation. It holds
    the source, media, instrument and probe systems, implements the core aspects of the photon
    packet life-cycle, and manages the iterative processes in the simulation (including phases,
    iterations and segments). Running a Monte Carlo simulation with SKIRT essentially comes down to
    constructing an instance of the MonteCarloSimulation class and invoking the setupAndRun()
    function on it. The runSimulation() function performs all segments of the simulation, calling
    functions such as runPrimaryEmission() and runSecondaryEmission() as required. These functions
    in turn invoke the performLifeCycle() function to trace photon packets through their complete
    life cycle, including emission, scattering events, peel-off towards the instruments, and
    registration of the contribution to the radiation field in each spatial cell crossed.

    The MonteCarloSimulation class also holds the non-discoverable \em config property, which is
    automatically set to an instance of the Configuration class. The setup() function of the config
    object is invoked at the very early stages of overall simulation setup, so that it can
    initialize its internal state to reflect the simulation configuration. As a result, it is safe
    for other simulation items to retrieve information from the config object during setup.

    <b>%Simulation phases, iterations and segments</b>

    In a first simulation phase, SKIRT determines the radiation field (RF) based on just the
    primary emission. In a second phase, it performs secondary emission using spectra calculated
    based on this RF.

    Some SKIRT media allow the density and/or the optical properties of the medium to depend on the
    local RF. These changes to the medium properties may sufficiently affect the RF to again cause
    a significant change in the medium properties, so that the RF must be calculated
    self-consistently by iterating over primary and/or secondary emission. We call this a \em
    dynamic \em medium \em state (DMS). A DMS update algorithm can be defined as part of the
    MaterialMix subclass that defines the complete behavior of a given material type, or it can be
    defined separately in a DynamicStateRecipe subclass, which can access any of the configured
    medium components. The latter option enables a single recipe to control different material
    types and it allows offering alternate recipes without complicating the material mix
    implementation.

    If the changes to the medium state significantly affect the medium opacity in the wavelength
    range of primary sources, iteration over primary emission is required. We call this a \em
    primary-dynamic \em medium \em state (PDMS). Radiative dust destruction is an example of a PDMS
    process.

    Alternatively, the media properties depending on the RF may be significant only for determining
    the secondary emission spectrum without a noticable effect during primary emission. As a
    result, there is no need for iteration over primary emission. We call this a \em
    secondary-dynamic \em medium \em state (SDMS). The hydrogen spin-flip transition is an example
    of an SDMS process because the opacity change at 21 cm does not affect primary emission at much
    shorter wavelengths (in this specific example, there is actually no need for iteration over the
    secondary emission either, but this might be different for other processes).

    Similarly, secondary emission may sufficiently affect the RF to cause a significant change in
    the secondary emission spectrum (which is calculated from that RF). If this is the case, the RF
    must be calculated self-consistently by iterating over the secondary emission. We call this \em
    dynamic \em secondary \em emission (DSE). One can consider DSE as a special case of SDMS where
    the change in the medium state is handled implicitly during the calculations rather than being
    explicitly stored. Dust emission in high-optical depth regions is an example of a DSE process.
    The dust temperature (or temperature probability distribution) is calculated on the fly while
    determining the emission spectrum from the RF.

    A simulation consecutively performs one or more \em segments as part of the relevant \em
    iterations and \em phases. A segment processes a set of photon packets, including emission from
    the relevant sources or media, absorption and scattering by the media, and updating the RF, the
    medium state, and/or the instruments as appropriate. We define the following segment types:

    - \f$\mathbf{P}\f$ -- primary emission: photon packets are launched from the configured
    (primary) source components.

    - \f$\mathbf{S}\f$ -- secondary emission: photon packets are launched from the configured
    medium components that support emission; this includes calculation of the secondary emission
    spectrum for each component based on the RF and the medium state in each spatial cell at the
    start of the segment.

    These segment identifiers are decorated by one or more of the following modifiers:

    - \f$p\f$: peel off a photon packet to the instruments for each emission and scattering event.

    - \f$r\f$: register the RF while tracing photon packets through the medium.

    - \f$(r)\f$: optionally register the RF if requested for probing (not needed for instruments).

    - \f$s\f$: update the PDMS based on the RF after processing all photon packets for the segment.

    - \f$(s)\f$: update the SDMS, if any, based on the RF after processing all photon packets for
    the segment.

    <b>%Simulation mode configuration</b>

    The user-configurable \em simulationMode enumeration sets the overall simulation mode. It
    determines the presence or absence of a transfer medium, the wavelength regime (oligochromatic
    or panchromatic), whether there is a secondary emission phase, and if so, which media types
    (dust and/or gas) are emitting photon packets. The Boolean user options \em
    iteratePrimaryEmission and \em iterateSecondaryEmission (plus \em includePrimaryEmission in the
    IterationOptions) further determine whether the simulation iterates over the radiation field to
    self-consistently calculate a DMS and/or DSE. The presence of DMS update algorithms is detected
    at run-time from the capabilities advertised by the configured media and recipes.

    An important aspect determined by the simulation mode is the simulation's wavelength regime,
    which can be oligochromatic or panchromatic. Oligochromatic simulations use just a few
    pre-defined, discrete wavelengths. They do no support kinematics (moving sources and/or media)
    because the wavelengths cannot shift away from the pre-defined values, and they do not support
    secondary emission by the transfer medium because the radiation field must be known across a
    wide spectrum to calculate the medium state and the resulting emission. Panchromatic
    simulations use a continuous range of wavelengths, lifting these limitations.

    The "extinction-only" simulation modes calculate the extinction of the primary radiation
    through the configured media, including the effects of absorption and scattering. There is no
    secondary emission, so these modes are meaningful only for wavelengths at which secondary
    sources (radiation from the media) can be neglected. The "emission" simulation modes (which
    require a panchromatic wavelength range) include secondary emission from dust and/or gas, in
    addition to the effects of absorption and scattering.

    The prefix \c Dust-, \c Gas-, or \c DustAndGas- in the emission simulation modes indicates
    which media types are emitting photon packets during secondary emission. Other media types can
    still be configured in the simulation and will contribute extinction during both primary and
    secondary emission.

    The simulation mode and iteration flags are set at the start of the configuration process.
    Their values have a significant impact on which options are allowed or required in the
    simulation's configuration.

    <b>Execution flow</b>

    The following table shows the flow of execution for each simulation mode. The arrows in the
    flow diagrams indicate forward progression (\f$\rightarrow\f$) and iteration
    (\f$\longleftarrow\f$).

    | | %Simulation mode | Flow diagram |
    |-|-----------------|--------------|
    | 1 | \c OligoNoMedium | \f$\mathbf{P}^p\f$ |
    | 2 | \c OligoExtinctionOnly | \f$\mathbf{P}^p_{(r)}\f$ |
    | 3 | \c LyaExtinctionOnly | \f$\mathbf{P}^p_{(r)}\f$ |
    | 4 | \c NoMedium | \f$\mathbf{P}^p\f$ |
    | 5 | \c ExtinctionOnly | \f$\mathbf{P}^p_{(r)}\f$ |
    | 6 | ... + \em iteratePrimaryEmission | \f$\overleftarrow{\mathbf{P}_{rs}} \;\rightarrow\; \mathbf{P}^p_{(r)}\f$ |
    | 7 | \c Dust-, \c Gas-, \c DustAndGasEmission | \f$\mathbf{P}^p_{r(s)} \;\rightarrow\; \mathbf{S}^p_{(r)}\f$ |
    | 8 | ... + \em iteratePrimaryEmission | \f$\overleftarrow{\mathbf{P}_{rs}} \;\rightarrow\; \mathbf{P}^p_{r(s)} \;\rightarrow\; \mathbf{S}^p_{(r)}\f$ |
    | 9 | ... + \em iterateSecondaryEmission | \f$\mathbf{P}^p_{r(s)} \;\rightarrow\; \overleftarrow{\mathbf{S}_{r(s)}} \;\rightarrow\; \mathbf{S}^p_{(r)}\f$ |
    | 10 | ... + \em iteratePrimaryEmission + \em iterateSecondaryEmission | \f$\overleftarrow{\mathbf{P}_{rs}} \;\rightarrow\; \textbf{}\mathbf{P}^p_{r(s)} \;\rightarrow\; \overleftarrow{\mathbf{S}_{r(s)}} \;\rightarrow\; \mathbf{S}^p_{(r)}\f$ |
    | 11 | ... + \em iterateSecondaryEmission + \em includePrimaryEmission | \f$\overleftarrow{\mathbf{P}_{r(s)} \;\rightarrow\; \mathbf{S}_{rs}} \;\rightarrow\; \mathbf{P}^p_{r(s)} \;\rightarrow\; \mathbf{S}^p_{(r)}\f$ |
    | 12 | ... + \em iteratePrimaryEmission + \em iterateSecondaryEmission + \em includePrimaryEmission | \f$\overleftarrow{\mathbf{P}_{rs}} \;\rightarrow\; \overleftarrow{\mathbf{P}_{r(s)} \;\rightarrow\; \mathbf{S}_{rs}} \;\rightarrow\; \mathbf{P}^p_{r(s)} \;\rightarrow\; \mathbf{S}^p_{(r)}\f$ |

    Rows 1-6 describe simulation modes that include primary emission only, without or with media.
    The last mode in this list (row 6) includes iteration over the PDMS during primary emission.
    Rows 7-12 describe modes that include dust and/or gas emission, possibly iterating over the RF
    during secondary emission to handle DSE (rows 9-12). If one or more of the media have a PDMS,
    there are modes supporting iteration over primary emission (rows 8, 10, 12) and/or over
    consecutive primary and secondary emission segments (rows 11, 12).

    If applicable, the simulation modes in rows 7-12 update the SDMS at the end of the relevant
    segments. It is meaningful to update an SDMS at the end of a segment that performs peel-off to
    the instruments because the update does not affect the outcome of a primary emission segment.

    <b>Photon life cycle</b>

    As mentioned above, each segment processes the life cycles for a set of photon packets,
    including emission and propagation through the medium and up to the instruments. The
    PhotonPacketOptions allow configuring four basic variations of the photon packet life cycle by
    enabling or disabling forced scattering and/or explicit absorption. Each of these variations
    comes with specific advantages or drawbacks, as follows:

    - With forced scattering. Forced scattering tends to reduce noise for simulations with low to
    limited optical depths, such as for most dust models on galaxy-wide scales. However, for each
    scattering event, it requires the calculation of the geometry and optical depth of the full
    path up to the model boundary, regardless of the location of the scattering event along the
    path.

    - Without forced scattering. In models with very intensive scattering, such as for Lyman-alpha
    line transfer, the photon cycle without forced scattering is often the better choice, because
    it avoids calculating the path geometry and optical depth beyond the scattering location.
    However, our implementation does not support storing the radiation field, which means this
    option cannot be used when the simulation includes secondary emission or dynamic state
    iteration.

    - Without explicit absorption. The default technique uses the extinction (sum of scattering and
    absorption) along a photon packet's path to locate the next interaction point. This requires
    the cumulative extinction optical depth to be a nondecreasing function of path length. It is
    thus not possible to handle negative extinction cross sections.

    - With explicit absorption. This technique instead used the scattering optical depth to locate
    the next interaction point. While the scattering cross section still must be nonnegative, this
    allows the extinction cross section to be negative. The latter is the case for materials that
    exhibit stimulated emission as soon as the magnitude of the negative absorption cross section
    is larger than the scattering cross section. As a drawback, this technique requires calculating
    both the scattering and absorption optical depths for the photon packet path.

    Refer to the performLifeCycle() function for more information on these life cycle variations.
    */
class MonteCarloSimulation : public Simulation
{
    /** The enumeration type indicating the simulation mode, which determines the overall structure
        of the simulation and its capabilities. The choice made for the simulation mode has a
        significant impact on which options are allowed or required in the simulation's
        configuration. See the class header for more information. */
    ENUM_DEF(SimulationMode, ExtinctionOnly)
        ENUM_VAL(SimulationMode, ExtinctionOnly, "Extinction only (no secondary emission)")
    ENUM_END()

    ITEM_CONCRETE(MonteCarloSimulation, Simulation, "a Monte Carlo simulation")

        PROPERTY_ENUM(simulationMode, SimulationMode, "the overall simulation mode")
        ATTRIBUTE_DEFAULT_VALUE(simulationMode, "ExtinctionOnly")
        ATTRIBUTE_INSERT(simulationMode, "simulationModeExtinctionOnly:Panchromatic,ExtinctionOnly;")

        PROPERTY_DOUBLE(numPackets, "the default number of photon packets launched per simulation segment")
        ATTRIBUTE_MIN_VALUE(numPackets, "[0")
        ATTRIBUTE_MAX_VALUE(numPackets, "1e19]")
        ATTRIBUTE_DEFAULT_VALUE(numPackets, "1e6")

        PROPERTY_ITEM(sourceSystem, SourceSystem, "the source system")
        ATTRIBUTE_DEFAULT_VALUE(sourceSystem, "SourceSystem")

        PROPERTY_ITEM(mediumSystem, MediumSystem, "the medium system")
        ATTRIBUTE_DEFAULT_VALUE(mediumSystem, "!NoMedium:MediumSystem;")
        ATTRIBUTE_REQUIRED_IF(mediumSystem, "!NoMedium")

        PROPERTY_ITEM(instrumentSystem, InstrumentSystem, "the instrument system")
        ATTRIBUTE_DEFAULT_VALUE(instrumentSystem, "InstrumentSystem")

        PROPERTY_ITEM(probeSystem, ProbeSystem, "the probe system")
        ATTRIBUTE_DEFAULT_VALUE(probeSystem, "ProbeSystem")

    ITEM_END()

    /** \fn numPackets
        The number of photon packets is specified as a double-precision floating point number
        rather than as a 64-bit integer to avoid implementing yet another discoverable property
        type. As a side benefit, one can use exponential notation to specify a large number of
        photon packets. Also, note that a double can exactly represent all integers up to 9e15. The
        maximum number of photon packets is somewhat arbitrarily set to 1e19 because that number is
        close to the maximum number representable with a 64-bit unsigned integer. */

    //============= Construction - Setup - Destruction =============

protected:
    /** This function performs setup for the complete simulation hierarchy. It calls the regular
        setup() function and notifies the probe system when setup has been completed. */
    void setupSimulation() override;

    /** This function constructs a SecondarySourceSystem object if the simulation configuration
        requires secondary emission. */
    void setupSelfAfter() override;

    //======== Getters for Non-Discoverable Properties =======

public:
    /** Returns the Configuration object for this simulation hierarchy. */
    Configuration* config() const;

    //======================== Other Functions =======================

protected:
    /** This function actually runs the simulation, assuming setup has been completed. It
        determines the appropriate execution flow from the simulation's configuration (see the
        table in the class header documentation) and invokes other functions in this class to
        consecutively perform all segments of the simulation, including primary source segments,
        secondary source segments, and any required iterations. For each segment, these invoked
        functions will tell the source system to prepare for launching photon packets (in serial
        code), and then cause photon packets to be launched and traced through their life cycles
        (in appropriately parallelized code depending on the run-time environment and the
        command-line options). */
    void runSimulation() override;

private:
    /** This function runs a final primary source emission segment including peel-off towards the
        instruments. It records radiation field contributions if the configuration requires it
        (e.g., because secondary emission must be calculated, or because the user configured probes
        to directly output radiation field information). If the configuration also includes
        emission, any secondary dynamic state updates are performed at the end of the segment.

        Using the notation of the table in the class header documentation, this function implements
        the execution flow \f$\mathbf{P}^p\f$, \f$\mathbf{P}^p_{(r)}\f$, or
        \f$\mathbf{P}^p_{r(s)}\f$, depending on the relevant configuration options. */
    void runPrimaryEmission();

    /** This function runs a final secondary source emission segment including peel-off towards the
        instruments. It records radiation field contributions if the configuration requires it
        (e.g., because the user configured probes to directly output radiation field information).

        Using the notation of the table in the class header documentation, this function implements
        the \f$\mathbf{S}^p_{(r)}\f$ execution flow. */
    void runSecondaryEmission();

    /** This function iteratively runs consecutive primary source emission segments to
        self-consistently calculate the radiation field taking into account any recipes for
        updating the primary dynamic medium state as a function of the radiation field.

        Each segment in the iteration tracks photon packets through the medium and subsequently
        updates the medium state based on the newly established radiation field. The number of
        photon packets launched for each segment can be configured as a multiplication factor on
        the default number of primary photon packets. The minimum and maximum number of iterations
        can also be specified as configuration options. Within these limits, the actual number of
        iterations performed is determined by convergence criteria defined by the configured
        dynamic medium state recipe(s).

        The photon packets emitted by this function do not perform peel-off to the instruments
        because the radiation field has not yet converged. After this function returns, a final
        segment of primary emission with peel-off must be performed by calling the
        runPrimaryEmission() function.

        Using the notation of the table in the class header documentation, this function implements
        the \f$\overleftarrow{\mathbf{P}_{rs}}\f$ execution flow. */
    void runPrimaryEmissionIterations();

    /** This function iteratively runs consecutive secondary source emission segments to
        self-consistently calculate the radiation field in models where the secondary emission
        spectrum and the total (primary plus secondary) radiation field significantly depend on
        each other. A typical example is a model where the dust medium is sufficiently opaque at
        infrared wavelengths to cause significant self-absorption and thus additional heating.

        Each segment in the iteration tracks photon packets through the medium to calculate the
        radiation field resulting from the secondary sources in the simulation. The subsequent
        segment then recalculates the secondary emission spectrum based on the newly established
        total (primary and secondary) radiation field. The number of photon packets launched for
        each segment can be configured as a multiplication factor on the default number of primary
        photon packets. The minimum and maximum number of iterations can also be specified as
        configuration options. Within these limits, the actual number of iterations performed is
        determined by user-configured convergence criteria. For dust media, convergence is reached
        when (a) the absorbed secondary luminosity is less than a given fraction of the absorbed
        primary luminosity, \em OR (b) the absorbed secondary luminosity has changed by less than a
        given fraction compared to the previous iteration.

        Any secondary dynamic state update algorithms in the simulation are performed at the end of
        each segment. Convergence is determined by the specific mechanism includes in each of these
        algorithms.

        The photon packets emitted by this function do not perform peel-off to the instruments
        because the radiation field has not yet converged. After this function returns, a final
        segment of secondary emission with peel-off must be performed by calling the
        runSecondaryEmission() function.

        Using the notation of the table in the class header documentation, this function implements
        the \f$\overleftarrow{\mathbf{S}_{r(s)}}\f$ execution flow. */
    void runSecondaryEmissionIterations();

    /** This function iteratively runs consecutive primary and secondary source emission segments
        to self-consistently calculate the radiation field in models that have both a dynamic
        medium state (the extinction properties depend on the radiation field) and dynamic
        secondary emission (the emission spectrum depends on the radiation field), and the user
        requested to include the primary emission in secondary emission iterations. This would be
        meaningful in case the changes in the radiation field caused by secondary emission could
        sufficiently influence the medium state to significantly affect the results of primary
        emission.

        The iteration is considered to be converged only when both the criteria for primary
        emission iteration and those for secondary emission iteration are satisfied; for more
        information see the runPrimaryEmissionIterations() and runSecondaryEmissionIterations()
        functions.

        Any secondary dynamic state update algorithms in the simulation are performed at the end of
        each primary segment, and the primary dynamic state update algorithms are performed at the
        end of each secondary emission segment. Convergence is determined by the specific mechanism
        includes in each of these algorithms.

        The photon packets emitted by this function do not perform peel-off to the instruments
        because the radiation field has not yet converged. After this function returns, a final
        segment of secondary emission with peel-off must be performed by calling the
        runSecondaryEmission() function.

        Using the notation of the table in the class header documentation, this function implements
        the \f$\overleftarrow{\mathbf{P}_{r(s)} \;\rightarrow\; \mathbf{S}_{rs}}\f$ execution flow.
        */
    void runMergedEmissionIterations();

    /** This function initializes the progress counter used in logprogress() for the specified
        segment and logs the number of photon packets to be processed. */
    void initProgress(string segment, size_t numTotal);

    /** This function logs a progress message for the segment specified in the initprogress()
        function if the previous message was issued at least some time interval ago. The function
        must be called regularly while processing photon packets. The argument specifies the number
        of photon packets processed. */
    void logProgress(size_t numDone);

    //======================== Data Members ========================

private:
    // non-discoverable simulation items
    Configuration* _config{new Configuration(this)};
    SecondarySourceSystem* _secondarySourceSystem{nullptr};  // constructed only when there is secondary emission

    // data members used by the XXXprogress() functions in this class
    string _segment;  // a string identifying the photon shooting segment for use in the log message
};

////////////////////////////////////////////////////////////////////

#endif
