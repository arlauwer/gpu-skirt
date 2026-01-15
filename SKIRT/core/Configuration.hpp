/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

#include "Range.hpp"
#include "SimulationItem.hpp"
class DisjointWavelengthGrid;
class SpatialCellLibrary;
class WavelengthDistribution;
class WavelengthGrid;

////////////////////////////////////////////////////////////////////

/** Configuration is a helper class that serves as a central clearing house for overall simulation
    configuration options including the simulation mode.

    Each MonteCarloSimulation holds a single Configuration object. During setup, it retrieves many
    properties and options from the simulation hierarchy, verifying consistency of the
    configuration and flagging any conflicts while doing so. Once this process has completed, the
    Configuration object offers getters for these retrieved properties to any of the other
    simulation items in the hierarchy. The setup() function of the Configuration object is invoked
    at the very early stages of simulation setup, so that it is safe for other simulation items to
    retrieve information from the Configuration object during setup.

    The Configuration class is based on SimulationItem so that it can be part of a simulation item
    hierarchy, however it is not discoverable because it is not intended to be selected or
    configured by the user. */
class Configuration : public SimulationItem
{
    //============= Construction - Setup - Destruction =============

public:
    /** This constructor creates a Configuration object that is hooked up as a child to the
        specified parent in the simulation hierarchy, so that it will automatically be deleted. The
        setup() function is \em not called by this constructor. */
    explicit Configuration(SimulationItem* parent);

protected:
    /** This function retrieves properties and options from the simulation hierarchy and stores the
        resulting values internally so that they can be returned by any of the getters with minimal
        overhead. During this process, the function also verifies the consistency of the simulation
        configuration, for example checking the configuration against the requirements of the
        selected simulation mode. If any conflicts are found, the function throws a fatal error. */
    void setupSelfBefore() override;

    /** This function logs some aspects of the configuration as information to the user. */
    void setupSelfAfter() override;

    //=========== Getters for configuration properties ============

public:
    // ----> wavelengths

    /** Returns true if the wavelength regime of the simulation is oligochromatic. */
    bool oligochromatic() const { return _oligochromatic; }

    /** Returns the total wavelength range of the primary sources in the simulation. For
        panchromatic simulations, this range is configured by the user in the source system. For
        oligochromatic simulations, the range includes the discrete source wavelengths used in the
        simulation, which are also user-configured in the source system. */
    Range sourceWavelengthRange() const { return _sourceWavelengthRange; }

    /** Returns a wavelength range that covers all wavelengths possibly used in the simulation for
        photon transport or for otherwise probing material properties (e.g. optical depth). This
        range includes the primary and secondary source wavelength ranges extended on both sides to
        accommodate a redshift or blueshift caused by kinematics corresponding to \f$v/c=1/3\f$. It
        also includes the range of the instrument wavelength grids and the wavelengths used for
        material normalization and material property probes. */
    Range simulationWavelengthRange() const;

    /** Returns a list of wavelengths that are explicitly or indirectly mentioned by the simulation
        configuration. This includes the characteristic wavelengths of all configured wavelength
        grids (for instruments, probes, radiation field or dust emission) and specific wavelengths
        used for normalization or probing. */
    vector<double> simulationWavelengths() const;

    /** Returns the wavelength grid to be used for an instrument or probe, given the wavelength
        grid configured locally for the calling instrument or probe (which may the null pointer to
        indicate that no local grid was configured). For oligochromatic simulations, the function
        always returns a wavelength grid with disjoint bins centered around the discrete source
        wavelengths used in the simulation. For panchromatic simulations, the function returns the
        provided local wavelength grid if it is non-null, and otherwise it returns the default
        instrument wavelength grid obtained from the instrument system. If both the provided local
        wavelength grid and the default instrument wavelength grid are the null pointer, the
        function throws a fatal error. */
    WavelengthGrid* wavelengthGrid(WavelengthGrid* localWavelengthGrid) const;

    /** For oligochromatic simulations, this function returns the wavelength bias distribution to
        be used by all primary sources. For panchromatic simulations, the function returns the null
        pointer. */
    WavelengthDistribution* oligoWavelengthBiasDistribution() { return _oligoWavelengthBiasDistribution; }

    // ----> probes

    /** Returns true if one of the Snapshot::getEntities() functions may be called for any of the
        snapshots associated with the imported sources and media in the simulation, implying that
        the snapshot must prebuild the required search data structures. In the current
        implementation, this happens only if the simulation includes one or more input model
        probes, i.e. instances of an InputModelProbe subclass. */
    bool snapshotsNeedGetEntities() const { return _snapshotsNeedGetEntities; }

    // ----> media

    /** Returns true if there is at least one medium component in the simulation, and false
        otherwise. */
    bool hasMedium() const { return _hasMedium; }

    /** Returns true if the Medium::generatePosition() function may be called for the media in the
        simulation. In the current implementation, this happens only if the simulation uses a
        VoronoiMeshSpatialGrid instance to discretize the spatial domain. If there are no media or
        the Medium::generatePosition() will never be called during this simulation, this function
        returns false. */
    bool mediaNeedGeneratePosition() const { return _mediaNeedGeneratePosition; }

    /** Returns true if the perceived photon packet wavelength equals the intrinsic photon packet
        wavelength for all spatial cells along the path of the packet. The following conditions
        cause this function to return false: Hubble expansion is enabled or some media may have a
        non-zero velocity in some cells. */
    bool hasConstantPerceivedWavelength() const { return _hasConstantPerceivedWavelength; }

    /** Returns true if the simulation has a exactly one medium component and the absorption and
        scattering cross sections for a photon packet traversing that medium component are
        spatially constant, so that the opacity in each crossed cell can be calculated by
        multiplying this constant cross section by the number density in the cell. Otherwise the
        function returns false.

        The following conditions cause this function to return false: Hubble expansion is enabled,
        there is more than one medium component, the medium may have a non-zero velocity in some
        cells, the medium has a variable material mix; the cross sections for some material mixes
        depend on extra medium state variables such as temperature or fragment weight factors. */
    bool hasSingleConstantSectionMedium() const { return _hasSingleConstantSectionMedium; }

    /** Returns true if the simulation has two or more medium components and the absorption and
        scattering cross sections for a photon packet traversing those medium components are
        spatially constant, so that the opacity in each crossed cell can be calculated by
        multiplying these constant cross sections by the corresponding number densities in the
        cell. Otherwise the function returns false.

        The following conditions cause this function to return false: Hubble expansion is enabled,
        some media may have a non-zero velocity in some cells, so that the perceived wavelength
        changes between cells; some media have a variable material mix; the cross sections for some
        material mixes depend on extra medium state variables such as temperature or fragment
        weight factors. */
    bool hasMultipleConstantSectionMedia() const { return _hasMultipleConstantSectionMedia; }

    /** Returns true if a scattering interaction for one or more media may emulate secondary
        emission, and false otherwise. */
    bool scatteringEmulatesSecondaryEmission() const { return _scatteringEmulatesSecondaryEmission; }

    /** Returns true if a scattering interaction for one or more media may adjust the wavelength of
        the interacting photon packet or may emulate secondary emission, and false otherwise. */
    bool needIndividualPeelOff() const { return _needIndividualPeelOff; }

    /** Returns true if some of the media in the simulation represent spheroidal (i.e.
        non-spherical) particles and require the corresponding treatment of polarization for
        scattering, absorption and emission, or false otherwise. If this function returns true, the
        hasPolarization() and hasMagneticField() functions return true as well. */
    bool hasSpheroidalPolarization() const { return _hasSpheroidalPolarization; }

    // ----> media sampling

    /** Returns the number of random spatial samples for determining density (or mass). */
    int numDensitySamples() const { return _numDensitySamples; }

    /** Returns the number of random spatial samples for determining other properties. */
    int numPropertySamples() const { return _numPropertySamples; }

    // ----> phases, iterations, number of packets

    /** Returns the number of photon packets launched per regular primary emission simulation
        segment. */
    double numPrimaryPackets() const { return _numPrimaryPackets; }

    // ----> photon cycle

    /** Returns the minimum weight reduction factor before a photon packet is terminated. */
    double minWeightReduction() const { return _minWeightReduction; }

    /** Returns the minimum number of forced scattering events before a photon packet is
        terminated. */
    int minScattEvents() const { return _minScattEvents; }

    /** Returns the fraction of path lengths sampled from a linear rather than an exponential
        distribution. */
    double pathLengthBias() const { return _pathLengthBias; }

    // ----> radiation field

    /** Returns true if the radiation field must be stored during the photon cycle, and false
        otherwise. */
    bool hasRadiationField() const { return _hasRadiationField; }

    /** Returns the wavelength grid to be used for storing the radiation field, or the null pointer
        if hasRadiationField() returns false. */
    DisjointWavelengthGrid* radiationFieldWLG() const { return _radiationFieldWLG; }

    //======================== Data Members ========================

private:
    // wavelengths
    bool _oligochromatic{false};
    Range _sourceWavelengthRange;
    WavelengthGrid* _defaultWavelengthGrid{nullptr};
    WavelengthDistribution* _oligoWavelengthBiasDistribution{nullptr};

    // probes
    bool _snapshotsNeedGetEntities{false};

    // media
    bool _hasMedium{false};
    bool _mediaNeedGeneratePosition{false};
    bool _hasConstantPerceivedWavelength{false};
    bool _hasSingleConstantSectionMedium{false};
    bool _hasMultipleConstantSectionMedia{false};
    bool _hasScatteringDispersion{false};
    bool _scatteringEmulatesSecondaryEmission{false};
    bool _needIndividualPeelOff{false};
    bool _hasSpheroidalPolarization{false};

    // media sampling
    int _numDensitySamples{100};
    int _numPropertySamples{1};

    // phases, iterations, number of packets
    double _numPrimaryPackets{0.};

    // photon cycle
    double _minWeightReduction{1e4};
    int _minScattEvents{0};
    double _pathLengthBias{0.5};

    // radiation field
    bool _hasRadiationField{false};
    DisjointWavelengthGrid* _radiationFieldWLG{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
