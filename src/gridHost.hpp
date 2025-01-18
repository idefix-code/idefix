// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef GRIDHOST_HPP_
#define GRIDHOST_HPP_

#include <vector>
#include "idefix.hpp"
#include "grid.hpp"
#include "input.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////////
/// The GridHost class is designed to store the grid data of the FULL computational domain (i.e. of
/// all of the MPI processes running) on the HOST. It comes handy to define initial conditions on
/// the Host and for output routines. Typical usage is to instantiate a GridHost from an existing
/// Grid instance, sync it and the use the grid information in Host routines.
/////////////////////////////////////////////////////////////////////////////////////////////////

class GridHost {
 public:
  std::array<IdefixHostArray1D<real>,3> x;     ///< geometrical central points
  std::array<IdefixHostArray1D<real>,3> xr; ///< cell right interface
  std::array<IdefixHostArray1D<real>,3> xl; ///< cell left interface
  std::array<IdefixHostArray1D<real>,3> dx; ///< cell width

  std::array<real,3> xbeg;                   ///< Beginning of grid
  std::array<real,3> xend;                   ///< End of grid

  std::array<int,3> np_tot;                  ///< total number of grid points
  std::array<int,3> np_int;                  ///< internal number of grid points

  std::array<int,3> nghost;                  ///< number of ghost cells
  std::array<BoundaryType,3> lbound;         ///< Boundary condition to the left
  std::array<BoundaryType,3> rbound;         ///< Boundary condition to the right

  bool haveAxis=false;    ///< Do we require a special treatment of the axis in spherical coords?
  bool isRegularCartesian; ///< whether the grid is regular and cartesian or not

  explicit GridHost(Grid&);   ///< Constructor from a corresponding Grid on the Device.
                              ///< (NB: this constructor does not sync any data)
  GridHost() = default;       ///< default constructor, should not be used explicitely.

  void MakeGrid(Input &);      ///< create grid coordinates from the input data.

  void SyncFromDevice();      ///< Synchronize this to the device Grid
  void SyncToDevice();      ///< Synchronize this from the device Grid

 private:
  Grid *grid;
};

#endif // GRIDHOST_HPP_
