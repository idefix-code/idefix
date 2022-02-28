// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef GRIDHOST_HPP_
#define GRIDHOST_HPP_

#include "idefix.hpp"
#include "grid.hpp"
#include "input.hpp"


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

  bool haveAxis=false;    // Do we require a special treatment of the axis in spherical coords?

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
