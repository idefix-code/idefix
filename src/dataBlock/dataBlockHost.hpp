// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef DATABLOCK_DATABLOCKHOST_HPP_
#define DATABLOCK_DATABLOCKHOST_HPP_

#include <vector>

#include "idefix.hpp"
#include "dataBlock.hpp"

class DataBlockHost {
 public:
  // Local grid information
  std::vector<IdefixArray1D<real>::HostMirror> x;   ///> geometrical central points
  std::vector<IdefixArray1D<real>::HostMirror> xr;  ///> cell right interface
  std::vector<IdefixArray1D<real>::HostMirror> xl;  ///> cell left interface
  std::vector<IdefixArray1D<real>::HostMirror> dx;  ///> cell width

  IdefixArray3D<real>::HostMirror dV;     ///> cell volume
  std::vector<IdefixArray3D<real>::HostMirror> A;   ///> cell right interface area

  IdefixArray4D<real>::HostMirror Vc;     ///> Main cell-centered primitive variables index
  #if MHD == YES
  IdefixArray4D<real>::HostMirror Vs;     ///> Main face-centered primitive variables index
  IdefixArray4D<real>::HostMirror Ve;     ///> Main edge-centered primitive variables index
  IdefixArray4D<real>::HostMirror J;      ///> Current (only when haveCurrent is enabled)

  IdefixArray3D<real>::HostMirror Ex1;    ///> x1 electric field
  IdefixArray3D<real>::HostMirror Ex2;    ///> x2 electric field
  IdefixArray3D<real>::HostMirror Ex3;    ///> x3 electric field

  #endif
  IdefixArray4D<real>::HostMirror Uc;     ///> Main cell-centered conservative variables
  IdefixArray3D<real>::HostMirror InvDt;


  std::vector<real> xbeg;                        ///> Beginning of dataBlock
  std::vector<real> xend;                        ///> End of dataBlock

  std::vector<int> np_tot;                       ///> total number of grid points
  std::vector<int> np_int;                       ///> internal number of grid points

  std::vector<int> nghost;                       ///> number of ghost cells

  std::vector<BoundaryType> lbound;           ///> Boundary condition to the left
  std::vector<BoundaryType> rbound;           ///> Boundary condition to the right

  std::vector<int> beg;                       ///> Begining of internal indices
  std::vector<int> end;                       ///> End of internal indices

  std::vector<int> gbeg;                      ///> Begining of local block in the grid (internal)
  std::vector<int> gend;                      ///> End of local block in the grid (internal)

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
