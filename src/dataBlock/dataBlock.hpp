// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef DATABLOCK_DATABLOCK_HPP_
#define DATABLOCK_DATABLOCK_HPP_

#include <vector>
#include <string>
#include "idefix.hpp"
#include "grid.hpp"
#include "gridHost.hpp"
#include "hydro.hpp"

#include "rkl.hpp"

#ifdef WITH_MPI
#include "mpi.hpp"
#endif

//TODO(lesurg) What is this standing for?
#define BOUNDARY_

class DataBlock {
 public:
  // Local grid information
  IdefixArray1D<real> x[3];    // geometrical central points
  IdefixArray1D<real> xr[3];   // cell right interface
  IdefixArray1D<real> xl[3];   // cell left interface
  IdefixArray1D<real> dx[3];   // cell width
  IdefixArray1D<real> xgc[3];  // cell geometrical cell center
  IdefixArray1D<real> rt;      // In spherical coordinates, gives \tilde{r}
  IdefixArray1D<real> sinx2m;      // In spherical coordinates,
                               // gives sin(th) at a j-1/2 interface
  IdefixArray1D<real> tanx2m;      // In spherical coordinates,
                               // gives tan(th) at a j-1/2 interface
  IdefixArray1D<real> sinx2;   // In spherical coordinates, gives sin(th) at the cell center
  IdefixArray1D<real> tanx2;   // In spherical coordinates, gives tan(th) at the cell center
  IdefixArray1D<real> dmu;     // In spherical coordinates,
                               // gives the \theta volume = fabs(cos(th_m) - cos(th_p))

  real xbeg[3];                // Beginning of datablocl
  real xend[3];                // End of datablock

  IdefixArray3D<real> dV;      // cell volume
  IdefixArray3D<real> A[3];    // cell right interface area

  int np_tot[3];               // total number of grid points
  int np_int[3];               // internal number of grid points

  int nghost[3];               // number of ghost cells
  BoundaryType lbound[3];      // Boundary condition to the left
  BoundaryType rbound[3];      // Boundary condition to the right

  bool haveAxis = false;       // This DataBlock contains points on the axis and a special treatment
                               // has been required for these.

  int beg[3];                  // Begining of internal indices
  int end[3];                  // End of internal indices

  int gbeg[3];                 // Begining of local block in the grid (internal)
  int gend[3];                 // End of local block in the grid (internal)

  real dt;                     // Current timestep
  real t;                      // Current time

  Grid *mygrid;

  #ifdef WITH_MPI
  Mpi mpi;                     // Mpi object when WITH_MPI is set
  #endif

  // The Hydro object attached to this datablock
  Hydro hydro;

  // The RKL object attached to this datablock
  RKLegendre rkl;

  // init from a Grid object
  void InitFromGrid(Grid &, Input &);

  void MakeGeometry();

  // Dump current datablock to a file for inspection
  void DumpToFile(std::string);

  // Return the number of cells who have Nans
  int CheckNan();

  // Evolve this DataBlock by dt
  void EvolveStage();

  // Reset the variables needed at each major integration Stage
  void ResetStage();

  DataBlock();

 private:
  void WriteVariable(FILE* , int , int *, char *, void*);
};

#endif // DATABLOCK_DATABLOCK_HPP_
