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

//////////////////////////////////////////////////////////////////////////////////////////////////
/// The DataBlockHost class is designed to store most of the information coming from an associated
/// DataBlock on the Host. It comes handy to define initial conditions on the Host and for output
/// routines. Similarly to the DataBlock, all of the arrays defined here contains information
/// attached to the local MPI sub-domain only. In particular grid-related of the DataBlockHost
/// only contains information of the local portion of the grid that belongs to the current MPI
/// process. The full grid is defined in Grid and GridHost classes.
/////////////////////////////////////////////////////////////////////////////////////////////////

class DataBlockHost {
 public:
  // Local grid information
  std::vector<IdefixArray1D<real>::HostMirror> x;   ///< geometrical central points
  std::vector<IdefixArray1D<real>::HostMirror> xr;  ///< cell right interface
  std::vector<IdefixArray1D<real>::HostMirror> xl;  ///< cell left interface
  std::vector<IdefixArray1D<real>::HostMirror> dx;  ///< cell width

  IdefixArray3D<real>::HostMirror dV;     ///< cell volume
  std::vector<IdefixArray3D<real>::HostMirror> A;   ///< cell right interface area

  IdefixArray4D<real>::HostMirror Vc;     ///< Main cell-centered primitive variables index
  #if MHD == YES
  IdefixArray4D<real>::HostMirror Vs;     ///< Main face-centered primitive variables index
  IdefixArray4D<real>::HostMirror Ve;     ///< Main edge-centered primitive variables index
  IdefixArray4D<real>::HostMirror J;      ///< Current (only when haveCurrent is enabled)

  IdefixArray3D<real>::HostMirror Ex1;    ///< x1 electric field
  IdefixArray3D<real>::HostMirror Ex2;    ///< x2 electric field
  IdefixArray3D<real>::HostMirror Ex3;    ///< x3 electric field

  #endif
  IdefixArray4D<real>::HostMirror Uc;     ///< Main cell-centered conservative variables
  IdefixArray3D<real>::HostMirror InvDt;  ///< Inverse of maximum timestep in each cell

  std::vector<IdefixArray2D<int>::HostMirror> coarseningLevel; ///< Grid coarsening level
                                                     ///< (only defined when coarsening
                                                     ///< is enabled)

  std::vector<bool> coarseningDirection;         ///< whether a coarsening is used in each direction


  std::vector<real> xbeg;                        ///> Beginning of dataBlock
  std::vector<real> xend;                        ///> End of dataBlock

  std::vector<int> np_tot;                  ///< total number of grid points
  std::vector<int> np_int;                  ///< internal number of grid points

  std::vector<int> nghost;                  ///< number of ghost cells

  std::vector<BoundaryType> lbound;           ///< Boundary condition to the left
  std::vector<BoundaryType> rbound;           ///< Boundary condition to the right

  std::vector<int> beg;                       ///< Begining of internal indices
  std::vector<int> end;                       ///< End of internal indices

  std::vector<int> gbeg;                      ///< Begining of local block in the grid (internal)
  std::vector<int> gend;                      ///< End of local block in the grid (internal)

  explicit DataBlockHost(DataBlock &);        ///< Constructor from a device datablock
                                              ///< (NB: does not sync any data)

  DataBlockHost() = default;                  ///< Default constructor

  void MakeVsFromAmag(IdefixHostArray4D<real> &); ///< Compute a face-centered mag. field in Vs from
                                                  ///< potential vector in argument

  void SyncToDevice();                            ///< Synchronize this to the device datablock
  void SyncFromDevice();                          ///< Synchronize this from the device datablock

  bool haveCurrent;                               ///< Whether the electrical current J is defined

  GridCoarsening haveGridCoarsening{GridCoarsening::disabled}; ///< Is grid coarsening enabled?

 private:
  // Data object to which we are the mirror
  DataBlock *data;
};

#endif // DATABLOCK_DATABLOCKHOST_HPP_
