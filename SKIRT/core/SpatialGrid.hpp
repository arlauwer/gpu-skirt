/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPATIALGRID_HPP
#define SPATIALGRID_HPP

#include "Box.hpp"
#include "Position.hpp"
#include "SimulationItem.hpp"
class PathSegmentGenerator;
class Random;
class SpatialGridPlotFile;

//////////////////////////////////////////////////////////////////////

/** The SpatialGrid class is an abstract base class for grids that tessellate the spatial domain of
    the simulation. Each position in the computational domain corresponds to a single spatial cell.
    A SpatialGrid subclass instance represents only purely geometric properties, i.e.\ it contains
    no information on the actual distribution of material over the grid. */
class SpatialGrid : public SimulationItem
{
    ITEM_ABSTRACT(SpatialGrid, SimulationItem, "a spatial grid")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function caches the simulation's random generator for use by subclasses. */
    void setupSelfBefore() override;

    //================ Functions that must be implemented in subclasses ===============

public:
    /** This function returns the number of cells in the grid. */
    virtual int numCells() const = 0;

    /** This function returns the bounding box that encloses the grid. */
    virtual Box boundingBox() const = 0;

    /** This function returns the volume of the cell with index \f$m\f$. */
    virtual double volume(int m) const = 0;

    /** This function returns the actual or approximate diagonal of the cell with index \f$m\f$.
        For cuboidal cells, the function returns the actual diagonal. For other geometric forms, it
        returns some approximate diagonal. */
    virtual double diagonal(int m) const = 0;

    /** This function returns the index \f$m\f$ of the cell that contains the position
        \f${\bf{r}}\f$. */
    virtual int cellIndex(Position bfr) const = 0;

    /** This function returns the central location of the cell with index \f$m\f$. */
    virtual Position centralPositionInCell(int m) const = 0;

    /** This function returns a random location from the cell with index \f$m\f$. */
    virtual Position randomPositionInCell(int m) const = 0;

    //======================== Other Functions =======================

protected:
    /** This function returns the simulation's random generator as a service to subclasses. */
    Random* random() const { return _random; }

    //======================== Data Members ========================

private:
    // data member initialized during setup
    Random* _random{nullptr};
};

//////////////////////////////////////////////////////////////////////

#endif
