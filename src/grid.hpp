#ifndef GRID_HPP
#define GRID_HPP
#include "idefix.hpp"

class Grid {
public:
    IdefixArray1D<real> x[3];      // geometrical central points
    IdefixArray1D<real>::HostMirror xH[3];
    IdefixArray1D<real> xr[3];     // cell right interface
    IdefixArray1D<real>::HostMirror xrH[3];
    IdefixArray1D<real> xl[3];     // cell left interface
    IdefixArray1D<real>::HostMirror xlH[3];
    IdefixArray1D<real> dx[3];     // cell width
    IdefixArray1D<real>::HostMirror dxH[3];

    IdefixArray3D<real> dV;     // cell volume
    IdefixArray3D<real>::HostMirror dVH;
    IdefixArray3D<real> A[3];      // cell right interface area
    IdefixArray3D<real>::HostMirror AH[3];


    // Constructor
    Grid(Input&);
    
};

#endif