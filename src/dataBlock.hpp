#ifndef DATABLOCK_HPP
#define DATABLOCK_HPP
#include "idefix.hpp"

class DataBlock {
public:
    // Local grid information
    IdefixArray1D<real> x[3];      // geometrical central points
    IdefixArray1D<real> xr[3];     // cell right interface
    IdefixArray1D<real> xl[3];     // cell left interface
    IdefixArray1D<real> dx[3];     // cell width

    IdefixArray3D<real> dV;     // cell volume
    IdefixArray3D<real> A[3];      // cell right interface area

    IdefixArray4D<real> Vc;     // Main cell-centered primitive variables index
    IdefixArray4D<real> Uc;     // Main cell-centered conservative variables
   
    int np_tot[3];                  // total number of grid points
    int np_int[3];                  // internal number of grid points

    int nghost[3];                  // number of ghost cells
    int lbound[3];                  // Boundary condition to the left
    int rbound[3];                  // Boundary condition to the right

    int beg[3];                     // Begining of internal indices
    int end[3];                     // End of internal indices

    int gbeg[3];                    // Begining of local block in the grid (internal)
    int gend[3];                    // End of local block in the grid (internal)

    // Constructor
    DataBlock(Grid &);

    // Copy constructor
    DataBlock(const DataBlock &);

    // Assignement operator
    DataBlock& operator=(const DataBlock&);


    DataBlock();
};



#endif