// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef GRID_HPP_
#define GRID_HPP_
#include "idefix.hpp"
#include "input.hpp"


class Grid {
 public:
  IdefixArray1D<real> x[3];    // geometrical central points
  IdefixArray1D<real> xr[3];   // cell right interface
  IdefixArray1D<real> xl[3];   // cell left interface
  IdefixArray1D<real> dx[3];   // cell width

  real xbeg[3];           // Beginning of grid
  real xend[3];           // End of grid

  int np_tot[3];          // total number of grid points (including ghosts)
  int np_int[3];          // internal number of grid points (excluding ghosts)

  int nghost[3];          // number of ghost cells
  BoundaryType lbound[3];          // Boundary condition to the left
  BoundaryType rbound[3];          // Boundary condition to the right

  bool haveAxis=false;    // Do we require a special treatment of the axis in spherical coords?

  // MPI data
  int nproc[3];           // Total number of procs in each direction
  int xproc[3];           // Coordinates of current proc in the array of procs

  #ifdef WITH_MPI
  MPI_Comm CartComm;
  #endif

  // Constructor
  explicit Grid(Input &);

  Grid();

 private:
  // Check if number is a power of 2
  bool isPow2(int);
  void makeDomainDecomposition();
};

#endif // GRID_HPP_
