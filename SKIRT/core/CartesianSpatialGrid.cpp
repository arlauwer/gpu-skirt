/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CartesianSpatialGrid.hpp"
#include "NR.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void CartesianSpatialGrid::setupSelfAfter()
{
    BoxSpatialGrid::setupSelfAfter();

    // initialize our local mesh arrays
    _Nx = _meshX->numBins();
    _Ny = _meshY->numBins();
    _Nz = _meshZ->numBins();
    _xv = _meshX->mesh() * (xmax() - xmin()) + xmin();
    _yv = _meshY->mesh() * (ymax() - ymin()) + ymin();
    _zv = _meshZ->mesh() * (zmax() - zmin()) + zmin();
}

//////////////////////////////////////////////////////////////////////

int CartesianSpatialGrid::numCells() const
{
    return _Nx * _Ny * _Nz;
}

//////////////////////////////////////////////////////////////////////

double CartesianSpatialGrid::volume(int m) const
{
    return box(m).volume();
}

//////////////////////////////////////////////////////////////////////

double CartesianSpatialGrid::diagonal(int m) const
{
    return box(m).diagonal();
}

//////////////////////////////////////////////////////////////////////

int CartesianSpatialGrid::cellIndex(Position bfr) const
{
    int i = NR::locateFail(_xv, bfr.x());
    int j = NR::locateFail(_yv, bfr.y());
    int k = NR::locateFail(_zv, bfr.z());
    if (i < 0 || j < 0 || k < 0)
        return -1;
    else
        return index(i, j, k);
}

//////////////////////////////////////////////////////////////////////

Position CartesianSpatialGrid::centralPositionInCell(int m) const
{
    return Position(box(m).center());
}

//////////////////////////////////////////////////////////////////////

Position CartesianSpatialGrid::randomPositionInCell(int m) const
{
    return random()->position(box(m));
}

//////////////////////////////////////////////////////////////////////

int CartesianSpatialGrid::index(int i, int j, int k) const
{
    return k + _Nz * j + _Nz * _Ny * i;
}

//////////////////////////////////////////////////////////////////////

Box CartesianSpatialGrid::box(int m) const
{
    int i = m / (_Nz * _Ny);
    int j = (m / _Nz) % _Ny;
    int k = m % _Nz;

    if (i < 0 || j < 0 || k < 0 || i >= _Nx || j >= _Ny || k >= _Nz)
        return Box();
    else
        return Box(_xv[i], _yv[j], _zv[k], _xv[i + 1], _yv[j + 1], _zv[k + 1]);
}

//////////////////////////////////////////////////////////////////////
