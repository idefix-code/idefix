// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef GRIDHOST_HPP_
#define GRIDHOST_HPP_

#include "idefix.hpp"


class GridHost {
 public:
  IdefixArray1D<real>::HostMirror x[3]; // geometrical central points
  IdefixArray1D<real>::HostMirror xr[3]; // cell right interface
  IdefixArray1D<real>::HostMirror xl[3]; // cell left interface
  IdefixArray1D<real>::HostMirror dx[3]; // cell width

  real xbeg[3];                   // Beginning of grid
  real xend[3];                   // End of grid

  int np_tot[3];                  // total number of grid points
  int np_int[3];                  // internal number of grid points

  int nghost[3];                  // number of ghost cells
  BoundaryType lbound[3];                  // Boundary condition to the left
  BoundaryType rbound[3];                  // Boundary condition to the right


  // Constructor
  explicit GridHost(Grid&);
  GridHost();

  // Actually make the grid
  void MakeGrid(Input &);

  // Sync from a device grid
  void SyncFromDevice();

  // Sync to a device grid
  void SyncToDevice();

 private:
  Grid *grid;
};

#endif // GRIDHOST_HPP_
