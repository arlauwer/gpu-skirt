/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEDIUMSYSTEM_HPP
#define MEDIUMSYSTEM_HPP

#include "MaterialMix.hpp"
#include "Medium.hpp"
#include "MediumState.hpp"
#include "PhotonPacketOptions.hpp"
#include "RadiationFieldOptions.hpp"
#include "SamplingOptions.hpp"
#include "SimulationItem.hpp"
#include "SpatialGrid.hpp"
#include "Table.hpp"
class Configuration;
class MaterialState;
class PhotonPacket;
class Random;
class ShortArray;
class SpatialGridPath;
class WavelengthGrid;

//////////////////////////////////////////////////////////////////////

/** An instance of the MediumSystem class represents a complete medium system, which is the
    superposition of one or more transfer media. Each individual medium represents a spatial
    density distribution and defines the material properties of the medium at each location. While
    the specific material properties may vary with location, the fundamental material type must be
    the same throughout the spatial domain for each medium.

    In addition to the media input model, the MediumSystem class includes the spatial grid that
    tessellates the spatial domain of the simulation into cells, and manages the medium state and
    the radiation field for each spatial cell in this grid. The class therefore plays a central
    role in the simulation and it offers quite a few different type of functions.

    <b>Overall medium configuration</b>

    These functions offer basic information on the spatial grid, the medium components, and the
    material types in the medium system.

    <b>%Medium state</b>

    The medium state for each spatial cell in the simulation includes a set of \em common state
    variables shared by all medium components and a set of \em specific state variables for each
    individual medium component. See the MediumState class for more information. The functions in
    this section allow access to the common state variables for a given spatial cell.

    <b>Low-level optical properties</b>

    These functions allow retrieving absorption, scattering and extinction cross sections for a
    given spatial cell and material type as a function of wavelength, assuming fixed, predefined
    values for any quantities other than wavelength (e.g., a default temperature, no polarization,
    no kinematics). The values returned by these low-level functions may be used only during setup
    and for probing.

    <b>High-level photon life cycle</b>

    These functions support the photon life cycle by offering a generic interface for operations
    including tracing photon paths through the medium and performing scattering events. The
    functions calculate opacities and other medium properties in a given spatial cell by passing
    the incoming photon packet and the full medium state to the appropriate material mix for each
    medium component. This allows proper treatment of polarization and kinematics and supports
    dependencies on temperature or custom special state variables.

    <b>Radiation field</b>

    These functions allow storing the radiation field during the photon life cycle and retrieving
    the results after a set of photon's have been processed. The contribution to the radation field
    for each spatial cell and for each wavelength in the simulation's radiation field wavelength
    grid is traced separately for primary and secondary sources. This avoids the need for repeating
    primary emission during dust-temperature convergence iterations. At all times, the sum of the
    primary and secondary contributions represents the radiation field to be used as input for
    calculations. There is a third, temporary table that serves as a target for storing the
    secondary radiation field so that the "stable" primary and secondary tables remain available
    for calculating secondary emission spectra while shooting secondary photons through the grid.

    <b>Indicative temperature</b>

    These functions determine an indicative temperature in a given spatial cell and for a given
    material type, averaged over medium components if applicable. Depending on the material type,
    the indicative temperature may be based on the radiation field or it may be derived from a
    value given in the input model. In any case, it does not reflect a physical quantity and it
    should be used only for setup and probing purposes.

    <b>Emission</b>

    These functions support the secondary source system by offering a generic interface for
    calculating the secondary emission properties in a given spatial cell and for a given material
    type, including luminosities, spectra and polarization. */
class MediumSystem : public SimulationItem
{
    ITEM_CONCRETE(MediumSystem, SimulationItem, "a medium system")
        ATTRIBUTE_TYPE_ALLOWED_IF(MediumSystem, "!NoMedium")

        PROPERTY_ITEM(photonPacketOptions, PhotonPacketOptions, "the photon packet options")
        ATTRIBUTE_DEFAULT_VALUE(photonPacketOptions, "PhotonPacketOptions")
        ATTRIBUTE_RELEVANT_IF(media, "!NoMedium")

        PROPERTY_ITEM(radiationFieldOptions, RadiationFieldOptions, "the radiation field options")
        ATTRIBUTE_DEFAULT_VALUE(radiationFieldOptions, "RadiationFieldOptions")
        ATTRIBUTE_RELEVANT_IF(radiationFieldOptions, "!NoMedium")

        PROPERTY_ITEM_LIST(media, Medium, "the transfer media")
        ATTRIBUTE_DEFAULT_VALUE(media, "GeometricMedium")
        ATTRIBUTE_REQUIRED_IF(media, "!NoMedium")

        PROPERTY_ITEM(samplingOptions, SamplingOptions, "the spatial grid sampling options")
        ATTRIBUTE_DEFAULT_VALUE(samplingOptions, "SamplingOptions")

        PROPERTY_ITEM(grid, SpatialGrid, "the spatial grid")
        ATTRIBUTE_DEFAULT_VALUE(grid,
                                "Dimension3:PolicyTreeSpatialGrid;Dimension2:Cylinder2DSpatialGrid;Sphere1DSpatialGrid")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates and stores initial state information for each spatial cell,
        including the cell volume and the number density for each medium as defined by the input
        model. If needed for the simulation's configuration, it also allocates one or two radiation
        field data tables that have a bin for each spatial cell in the simulation and for each bin
        in the wavelength grid returned by the Configuration::radiationFieldWLG() function. */
    void setupSelfAfter() override;

    //=============== Overall medium configuration ===================

public:
    /** This function returns the number of media in the medium system. The returned value is valid
        only after setup has been performed. */
    int numMedia() const;

    /** This function returns the number of cells in the spatial grid held by the medium system.
        The returned value is valid only after setup has been performed. */
    int numCells() const;

    /** This function returns the material mix corresponding to the medium component with index
        \f$h\f$ in spatial cell with index \f$m\f$. */
    const MaterialMix* mix(int m, int h) const;

    /** This function returns true if at least one of the media in the medium system has the
        specified fundamental material type (i.e. dust, electrons, or gas). */
    bool hasMaterialType(MaterialMix::MaterialType type) const;

    /** This function returns true if at least one of the media in the medium system contains dust.
        */
    bool hasDust() const { return hasMaterialType(MaterialMix::MaterialType::Dust); }

    /** This function returns true if at least one of the media in the medium system contains
        electrons. */
    bool hasElectrons() const { return hasMaterialType(MaterialMix::MaterialType::Electrons); }

    /** This function returns true if at least one of the media in the medium system contains gas.
        */
    bool hasGas() const { return hasMaterialType(MaterialMix::MaterialType::Gas); }

    /** This function returns true if the medium component with index \f$h\f$ has the specified
        fundamental material type (i.e. dust, electrons, or gas). */
    bool isMaterialType(MaterialMix::MaterialType type, int h) const;

    /** This function returns true if the medium component with index \f$h\f$ contains dust. */
    bool isDust(int h) const { return isMaterialType(MaterialMix::MaterialType::Dust, h); }

    /** This function returns true if the medium component with index \f$h\f$ contains electrons.
        */
    bool isElectrons(int h) const { return isMaterialType(MaterialMix::MaterialType::Electrons, h); }

    /** This function returns true if the medium component with index \f$h\f$ contains gas. */
    bool isGas(int h) const { return isMaterialType(MaterialMix::MaterialType::Gas, h); }

    /** This function returns a list of indices \f$h\f$ for media components that contain dust. */
    const vector<int>& dustMediumIndices() const { return _dust_hv; }

    /** This function returns a list of indices \f$h\f$ for media components that contain gas. */
    const vector<int>& gasMediumIndices() const { return _gas_hv; }

    /** This function returns a list of indices \f$h\f$ for media components that contain electrons. */
    const vector<int>& electronMediumIndices() const { return _elec_hv; }

    //=============== Input model ===================

public:
    /** This function returns the total mass density of all dust medium components at the specified
        position in the input model (i.e. directly querying the medium components rather than the
        gridded medium state). */
    double dustMassDensity(Position bfr) const;

    /** This function returns the total number density of all electron medium components at the specified
        position in the input model (i.e. directly querying the medium components rather than the
        gridded medium state). */
    double electronNumberDensity(Position bfr) const;

    /** This function returns the total number density of all gas medium components at the specified
        position in the input model (i.e. directly querying the medium components rather than the
        gridded medium state). */
    double gasNumberDensity(Position bfr) const;

    //=============== Medium state ===================

public:
    /** This function returns the volume of the spatial cell with index \f$m\f$. */
    double volume(int m) const;

    /** This function returns the total mass density of all medium components in spatial cell with
        index \f$m\f$. */
    double massDensity(int m) const;

    /** This function returns the total mass density of all dust medium components in spatial cell
        with index \f$m\f$. */
    double dustMassDensity(int m) const;

    /** This function returns the total number density of all electron medium components in spatial
        cell with index \f$m\f$. */
    double electronNumberDensity(int m) const;

    /** This function returns the total number density of all gas medium components in spatial cell
        with index \f$m\f$. */
    double gasNumberDensity(int m) const;

    /** This function returns the number density of the medium component with index \f$h\f$ in
        spatial cell with index \f$m\f$. */
    double numberDensity(int m, int h) const;

    /** This function returns the mass density of the medium component with index \f$h\f$ in
        spatial cell with index \f$m\f$. */
    double massDensity(int m, int h) const;

    /** This function returns the metallicity \f$Z\f$ of the medium component with index \f$h\f$ in
        spatial cell with index \f$m\f$. If the specified medium component does not have the
        metallicity specific state variable, the behavior of this function is undefined. */
    double metallicity(int m, int h) const;

    /** This function returns the temperature \f$T\f$ of the medium component with index \f$h\f$ in
        spatial cell with index \f$m\f$. If the specified medium component does not have the
        temperature specific state variable, the behavior of this function is undefined. */
    double temperature(int m, int h) const;

    /** This function returns the value of the custom specific state variable with index \f$i\f$ of
        the medium component with index \f$h\f$ in the spatial cell with index \f$m\f$. If the
        specified medium component does not have a custom variable with the specified index, the
        behavior of this function is undefined. */
    double custom(int m, int h, int i) const;

    //=============== Low-level optical properties ===================

public:
    /** This function returns the absorption opacity \f$k_h^\text{abs}\f$ at wavelength
        \f$\lambda\f$ of the medium component with index \f$h\f$ in spatial cell with index
        \f$m\f$. Because no photon packet is provided, default values are used for any relevant
        incoming photon packet properties. For example, the radiation is assumed to be unpolarized.
        */
    double opacityAbs(double lambda, int m, int h) const;

    /** This function returns the scattering opacity \f$k_h^\text{sca}\f$ at wavelength
        \f$\lambda\f$ of the medium component with index \f$h\f$ in spatial cell with index
        \f$m\f$. Because no photon packet is provided, default values are used for any relevant
        incoming photon packet properties. For example, the radiation is assumed to be unpolarized.
        */
    double opacitySca(double lambda, int m, int h) const;

    /** This function returns the extinction opacity \f$k_h^\text{ext}\f$ at wavelength
        \f$\lambda\f$ of the medium component with index \f$h\f$ in spatial cell with index
        \f$m\f$. Because no photon packet is provided, default values are used for any relevant
        incoming photon packet properties. For example, the radiation is assumed to be unpolarized.
        */
    double opacityExt(double lambda, int m, int h) const;

    /** This function returns the absorption opacity \f$k^\text{abs}=\sum_h k_h^\text{abs}\f$
        summed over all medium components with the specified material type at wavelength
        \f$\lambda\f$ in spatial cell with index \f$m\f$. Because no photon packet is provided,
        default values are used for any relevant incoming photon packet properties. For example,
        the radiation is assumed to be unpolarized. */
    double opacityAbs(double lambda, int m, MaterialMix::MaterialType type) const;

    /** This function returns the extinction opacity \f$k^\text{ext}=\sum_h k_h^\text{ext}\f$
        summed over all medium components with the specified material type at wavelength
        \f$\lambda\f$ in spatial cell with index \f$m\f$. Because no photon packet is provided,
        default values are used for any relevant incoming photon packet properties. For example,
        the radiation is assumed to be unpolarized. */
    double opacityExt(double lambda, int m, MaterialMix::MaterialType type) const;

    /** This function returns the extinction opacity \f$k^\text{ext}=\sum_h k_h^\text{ext}\f$
        summed over all medium components at wavelength \f$\lambda\f$ in spatial cell with index
        \f$m\f$. Because no photon packet is provided, default values are used for any relevant
        incoming photon packet properties. For example, the radiation is assumed to be unpolarized.
        */
    double opacityExt(double lambda, int m) const;

    //=============== High-level photon life cycle ===================

private:
    /** This function returns the absorption opacity \f$k^\text{abs}=\sum_h k_h^\text{abs}\f$
        summed over all medium components at wavelength \f$\lambda\f$ in spatial cell with index
        \f$m\f$, where applicable taking into account the properties of the specified incoming
        photon packet (for example, its polarization state). */
    double opacityAbs(double lambda, int m, const PhotonPacket* pp) const;

    /** This function returns the scattering opacity \f$k^\text{sca}=\sum_h k_h^\text{sca}\f$
        summed over all medium components at wavelength \f$\lambda\f$ in spatial cell with index
        \f$m\f$, where applicable taking into account the properties of the specified incoming
        photon packet (for example, its polarization state). */
    double opacitySca(double lambda, int m, const PhotonPacket* pp) const;

    /** This function returns the extinction opacity \f$k^\text{ext}=\sum_h k_h^\text{ext}\f$
        summed over all medium components at wavelength \f$\lambda\f$ in spatial cell with index
        \f$m\f$, where applicable taking into account the properties of the specified incoming
        photon packet (for example, its polarization state). */
    double opacityExt(double lambda, int m, const PhotonPacket* pp) const;

    //======================== Data Members ========================

private:
    Configuration* _config{nullptr};

    // relevant for any simulation mode that includes a medium
    int _numCells{0};  // index m
    int _numMedia{0};  // index h
    bool _mixPerCell{false};
    vector<const MaterialMix*> _mixv;  // material mixes; indexed on h, or on m and h if mixPerCell is true
    MediumState _state;                // state info for each cell and each medium

    // cached info relevant for any simulation mode that includes a medium
    vector<int> _dust_hv;  // a list of indices for media components containing dust
    vector<int> _gas_hv;   // a list of indices for media components containing gas
    vector<int> _elec_hv;  // a list of indices for media components containing electrons
    vector<int> _pdms_hv;  // a list of indices for media components with a primary dynamic medium state
    vector<int> _sdms_hv;  // a list of indices for media components with a secondary dynamic medium state

    // relevant for any simulation mode that stores the radiation field
    WavelengthGrid* _wavelengthGrid{0};  // index ell
    // each radiation field table has an entry for each cell and each wavelength (indexed on m,ell)
    // - the sum of rf1 and rf2 represents the stable radiation field to be used as input for regular calculations
    // - rf2c serves as a target for storing the secondary radiation field so that rf1+rf2 remain available for
    //   calculating secondary emission spectra while already shooting photons through the grid
    Table<2> _rf1;  // radiation field from primary sources

    // relevant for any simulation mode that includes dust emission
    int _numDustEmissionWavelengths{0};
};

////////////////////////////////////////////////////////////////

#endif
