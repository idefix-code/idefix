#ifndef GRIDHOST_HPP
#define GRIDHOST_HPP
#include "idefix.hpp"


class GridHost {
public:
    IdefixArray1D<real>::HostMirror x[3]; // geometrical central points
    IdefixArray1D<real>::HostMirror xr[3]; // cell right interface
    IdefixArray1D<real>::HostMirror xl[3]; // cell left interface
    IdefixArray1D<real>::HostMirror dx[3]; // cell width


    int np_tot[3];                  // total number of grid points
    int np_int[3];                  // internal number of grid points

    int nghost[3];                  // number of ghost cells
    BoundaryType lbound[3];                  // Boundary condition to the left
    BoundaryType rbound[3];                  // Boundary condition to the right

    
    // Constructor
    GridHost(Grid&);
    GridHost();

    // Actually make the grid
    void MakeGrid(Input &);

    // Sync from a device grid
    void SyncFromDevice();

    // Sync to a device grid
    void SyncToDevice();

private:
    Grid grid;
};

#endif