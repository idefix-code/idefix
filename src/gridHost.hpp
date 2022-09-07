// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
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
  std::vector<IdefixArray1D<real>::HostMirror> x;     ///< geometrical central points
  std::vector<IdefixArray1D<real>::HostMirror> xr; ///< cell right interface
  std::vector<IdefixArray1D<real>::HostMirror> xl; ///< cell left interface
  std::vector<IdefixArray1D<real>::HostMirror> dx; ///< cell width

  std::vector<real> xbeg;                   ///< Beginning of grid
  std::vector<real> xend;                   ///< End of grid

  std::vector<int> np_tot;                  ///< total number of grid points
  std::vector<int> np_int;                  ///< internal number of grid points

  std::vector<int> nghost;                  ///< number of ghost cells
  std::vector<BoundaryType> lbound;         ///< Boundary condition to the left
  std::vector<BoundaryType> rbound;         ///< Boundary condition to the right

  bool haveAxis=false;    ///< Do we require a special treatment of the axis in spherical coords?

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
