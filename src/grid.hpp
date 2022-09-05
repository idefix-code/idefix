// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef GRID_HPP_
#define GRID_HPP_
#include <vector>
#include "idefix.hpp"
#include "input.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////////
/// The Grid class is designed to store the grid data of the FULL computational domain on the Host
/// (i.e. of all of the MPI processes running).
/// The domain decomposition is performed by the child instances of the DataBlock class, which are
/// built on a Grid instance, using the cartesian MPI communicator of that Grid. Note that all of
/// the arrays of a Grid are on the device. If a Host access is needed,
/// it is recommended to use a GridHost instance from this Grid, and sync it.
/////////////////////////////////////////////////////////////////////////////////////////////////

class Grid {
 public:
  std::vector<IdefixArray1D<real>> x;    ///< geometrical central points
  std::vector<IdefixArray1D<real>> xr;   ///< cell right interface
  std::vector<IdefixArray1D<real>> xl;   ///< cell left interface
  std::vector<IdefixArray1D<real>> dx;   ///< cell width

  std::vector<real> xbeg;           ///< Beginning of grid
  std::vector<real> xend;           ///< End of grid

  std::vector<int> np_tot;          ///< total number of grid points (including ghosts)
  std::vector<int> np_int;          ///< internal number of grid points (excluding ghosts)

  std::vector<int> nghost;          ///< number of ghost cells
  std::vector<BoundaryType> lbound;          ///< Boundary condition to the left
  std::vector<BoundaryType> rbound;          ///< Boundary condition to the right

  bool haveAxis{false};    ///< Do we require a special treatment of the axis in spherical coords?


  GridCoarsening haveGridCoarsening{GridCoarsening::disabled}; ///< Is grid coarsening enabled?
  std::vector<bool> coarseningDirection;  ///< whether a coarsening is used in each direction

  // MPI data
  std::vector<int> nproc;           ///</< Total number of procs in each direction
  std::vector<int> xproc;           ///</< Coordinates of current proc in the array of procs

  #ifdef WITH_MPI
  MPI_Comm CartComm;                ///< Cartesian communicator for the planned domain decomposition
  MPI_Comm AxisComm;                ///< Cartesian communicator to exchange data accross the axis
                                    ///< (when applicable)
  #endif

  // Constructor
  explicit Grid(Input &);
  void ShowConfig();

  Grid() = default;

 private:
  // Check if number is a power of 2
  bool isPow2(int);
  void makeDomainDecomposition();
};

#endif // GRID_HPP_
