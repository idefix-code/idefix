// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
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
  std::array<IdefixHostArray1D<real>,3> x;   ///< geometrical central points
  std::array<IdefixHostArray1D<real>,3> xr;  ///< cell right interface
  std::array<IdefixHostArray1D<real>,3> xl;  ///< cell left interface
  std::array<IdefixHostArray1D<real>,3> dx;  ///< cell width

  IdefixHostArray3D<real> dV;     ///< cell volume
  std::array<IdefixHostArray3D<real>,3> A;   ///< cell right interface area

  IdefixHostArray4D<real> Vc;     ///< Main cell-centered primitive variables index

  bool haveDust{false};
  std::vector<IdefixHostArray4D<real>> dustVc; ///< Cell-centered primitive variables index for dust

  #if MHD == YES
  IdefixHostArray4D<real> Vs;     ///< Main face-centered primitive variables index
  IdefixHostArray4D<real> Ve;     ///< Main edge-centered primitive variables index
  IdefixHostArray4D<real> J;      ///< Current (only when haveCurrent is enabled)

  IdefixHostArray3D<real> Ex1;    ///< x1 electric field
  IdefixHostArray3D<real> Ex2;    ///< x2 electric field
  IdefixHostArray3D<real> Ex3;    ///< x3 electric field

  #endif
  IdefixHostArray4D<real> Uc;     ///< Main cell-centered conservative variables
  IdefixHostArray3D<real> InvDt;  ///< Inverse of maximum timestep in each cell

  std::array<IdefixHostArray2D<int>,3> coarseningLevel; ///< Grid coarsening level
                                                     ///< (only defined when coarsening
                                                     ///< is enabled)

  std::array<bool,3> coarseningDirection;       ///< whether a coarsening is used in each direction


  std::array<real,3> xbeg;                      ///> Beginning of dataBlock
  std::array<real,3> xend;                      ///> End of dataBlock

  std::array<int,3> np_tot;                  ///< total number of grid points
  std::array<int,3> np_int;                  ///< internal number of grid points

  std::array<int,3> nghost;                  ///< number of ghost cells

  std::array<BoundaryType,3> lbound;           ///< Boundary condition to the left
  std::array<BoundaryType,3> rbound;           ///< Boundary condition to the right

  std::array<int,3> beg;                       ///< Begining of internal indices
  std::array<int,3> end;                       ///< End of internal indices

  std::array<int,3> gbeg;                      ///< Begining of local block in the grid (internal)
  std::array<int,3> gend;                      ///< End of local block in the grid (internal)

  real dt;                     ///< Current timestep
  real t;                      ///< Current time

  explicit DataBlockHost(DataBlock &);        ///< Constructor from a device datablock
                                              ///< (NB: does not sync any data)

  DataBlockHost() = default;                  ///< Default constructor

  // The Planetary system (actually a copy from the dataBlock)
  bool haveplanetarySystem{false};
  PlanetarySystem* planetarySystem;


  void MakeVsFromAmag(IdefixHostArray4D<real> &); ///< Compute a face-centered mag. field in Vs from
                                                  ///< potential vector in argument

  void SyncToDevice();                            ///< Synchronize this to the device datablock
  void SyncFromDevice();                          ///< Synchronize this from the device datablock

  bool haveCurrent;                               ///< Whether the electrical current J is defined

  GridCoarsening haveGridCoarsening{GridCoarsening::disabled}; ///< Is grid coarsening enabled?

  // Data object to which we are the mirror
  DataBlock *data;
};

#endif // DATABLOCK_DATABLOCKHOST_HPP_
