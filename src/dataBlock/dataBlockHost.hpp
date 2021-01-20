// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef DATABLOCK_DATABLOCKHOST_HPP_
#define DATABLOCK_DATABLOCKHOST_HPP_

#include "idefix.hpp"
#include "dataBlock.hpp"

class DataBlockHost {
 public:
  // Local grid information
  IdefixArray1D<real>::HostMirror x[3];   // geometrical central points
  IdefixArray1D<real>::HostMirror xr[3];  // cell right interface
  IdefixArray1D<real>::HostMirror xl[3];  // cell left interface
  IdefixArray1D<real>::HostMirror dx[3];  // cell width

  IdefixArray3D<real>::HostMirror dV;     // cell volume
  IdefixArray3D<real>::HostMirror A[3];   // cell right interface area

  IdefixArray4D<real>::HostMirror Vc;     // Main cell-centered primitive variables index
  IdefixArray4D<real>::HostMirror Vs;     // Main face-centered primitive variables index
  IdefixArray4D<real>::HostMirror J;      // Current (only when haveCurrent is enabled)
  IdefixArray4D<real>::HostMirror Uc;     // Main cell-centered conservative variables
  IdefixArray3D<real>::HostMirror InvDt;

  IdefixArray3D<real>::HostMirror Ex1;    // x1 electric field
  IdefixArray3D<real>::HostMirror Ex2;    // x2 electric field
  IdefixArray3D<real>::HostMirror Ex3;    // x3 electric field


  real xbeg[3];                           // Beginning of dataBlock
  real xend[3];                           // End of dataBlock

  int np_tot[3];                          // total number of grid points
  int np_int[3];                          // internal number of grid points

  int nghost[3];                          // number of ghost cells

  BoundaryType lbound[3];                 // Boundary condition to the left
  BoundaryType rbound[3];                 // Boundary condition to the right

  int beg[3];                             // Begining of internal indices
  int end[3];                             // End of internal indices

  int gbeg[3];                            // Begining of local block in the grid (internal)
  int gend[3];                            // End of local block in the grid (internal)

  // Constructor
  explicit DataBlockHost(DataBlock &);

  // Default constructor
  DataBlockHost();

  // Construct a face-centered field from potential vector
  void MakeVsFromAmag(IdefixHostArray4D<real> &);

  // Synchronisation routines
  void SyncToDevice();
  void SyncFromDevice();

  // Whether or not the current is defined
  bool haveCurrent;


 private:
  // Data object to which we are the mirror
  DataBlock *data;
};

#endif // DATABLOCK_DATABLOCKHOST_HPP_
