#ifndef GRIDHOST_HPP
#define GRIDHOST_HPP
#include "idefix.hpp"


class GridHost {
public:
    IdefixArray1D<real>::HostMirror x[3]; // geometrical central points
    IdefixArray1D<real>::HostMirror xr[3]; // cell right interface
    IdefixArray1D<real>::HostMirror xl[3]; // cell left interface
    IdefixArray1D<real>::HostMirror dx[3]; // cell width

    IdefixArray3D<real>::HostMirror dV;     // cell volume
    IdefixArray3D<real>::HostMirror A[3];      // cell right interface area

    int np_tot[3];                  // total number of grid points
    int np_int[3];                  // internal number of grid points

    int nghost[3];                  // number of ghost cells
    int lbound[3];                  // Boundary condition to the left
    int rbound[3];                  // Boundary condition to the right

    int lbeg[3];                    // Begining of local array
    int lend[3];                    // end of local array
    
    // Constructor
    GridHost(Grid&, Input&);

    // Actually make the grid
    void MakeGrid(Grid &, Input &);

private:
    Grid grid;
};

#endif