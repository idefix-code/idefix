#ifndef DATABLOCK_HPP
#define DATABLOCK_HPP
#include "idefix.hpp"

#define     BOUNDARY_

#if MHD == YES
// Forward declaration of electromotive force
class ElectroMotiveForce;
#endif

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
    IdefixArray4D<real> Vs;     // Main face-centered varariables
    IdefixArray4D<real> Uc;     // Main cell-centered conservative variables

    // Required by time integrator
    IdefixArray4D<real> V0;
    IdefixArray3D<real> InvDtHyp;
    IdefixArray3D<real> InvDtPar;

    // Required by physics
    IdefixArray4D<real> PrimL;
    IdefixArray4D<real> PrimR;
    IdefixArray4D<real> FluxRiemann;
    

   
    int np_tot[3];                  // total number of grid points
    int np_int[3];                  // internal number of grid points

    int nghost[3];                  // number of ghost cells
    BoundaryType lbound[3];                  // Boundary condition to the left
    BoundaryType rbound[3];                  // Boundary condition to the right

    int beg[3];                     // Begining of internal indices
    int end[3];                     // End of internal indices

    int gbeg[3];                    // Begining of local block in the grid (internal)
    int gend[3];                    // End of local block in the grid (internal)

#if MHD == YES
    ElectroMotiveForce emf;
#endif
    // init from a Grid object
    void InitFromGrid(Grid &);

    // Copy constructor
    DataBlock(const DataBlock &);

    // Assignement operator
    DataBlock& operator=(const DataBlock&);


    DataBlock();


};

#if MHD == YES
#include "electromotiveforce.hpp"
#endif

#endif