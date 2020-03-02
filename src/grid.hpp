#ifndef GRID_HPP
#define GRID_HPP
#include "idefix.hpp"

class Grid {
public:
    IdefixArray1D<real> x[3];      // geometrical central points
    IdefixArray1D<real> xr[3];     // cell right interface
    IdefixArray1D<real> xl[3];     // cell left interface
    IdefixArray1D<real> dx[3];     // cell width


    int np_tot[3];                  // total number of grid points
    int np_int[3];                  // internal number of grid points

    int nghost[3];                  // number of ghost cells
    int lbound[3];                  // Boundary condition to the left
    int rbound[3];                  // Boundary condition to the right

    // Constructor
    Grid(Input &);

    Grid();



};

#endif