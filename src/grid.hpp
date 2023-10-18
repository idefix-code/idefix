// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef GRID_HPP_
#define GRID_HPP_
#include <vector>
#include <memory>
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
class SubGrid;    // Forward class declaration

class Grid {
 public:
  std::array<IdefixArray1D<real>,3> x;    ///< geometrical central points
  std::array<IdefixArray1D<real>,3> xr;   ///< cell right interface
  std::array<IdefixArray1D<real>,3> xl;   ///< cell left interface
  std::array<IdefixArray1D<real>,3> dx;   ///< cell width

  std::array<real,3> xbeg;           ///< Beginning of grid
  std::array<real,3> xend;           ///< End of grid

  std::array<int,3> np_tot;          ///< total number of grid points (including ghosts)
  std::array<int,3> np_int;          ///< internal number of grid points (excluding ghosts)

  std::array<int,3> nghost;          ///< number of ghost cells
  std::array<BoundaryType,3> lbound;          ///< Boundary condition to the left
  std::array<BoundaryType,3> rbound;          ///< Boundary condition to the right

  bool haveAxis{false};    ///< Do we require a special treatment of the axis in spherical coords?
  bool isRegularCartesian{true}; ///< whether the grid is regular and cartesian (true) or not

  GridCoarsening haveGridCoarsening{GridCoarsening::disabled}; ///< Is grid coarsening enabled?
  std::array<bool,3> coarseningDirection;  ///< whether a coarsening is used in each direction

  // MPI data
  std::array<int,3> nproc;           ///</< Total number of procs in each direction
  std::array<int,3> xproc;           ///</< Coordinates of current proc in the array of procs

  #ifdef WITH_MPI
  MPI_Comm CartComm;                ///< Cartesian communicator for the planned domain decomposition
  MPI_Comm AxisComm;                ///< Cartesian communicator to exchange data accross the axis
                                    ///< (when applicable)
  #endif

  // Constructor
  explicit Grid(Input &);
  explicit Grid(SubGrid *);

  void ShowConfig();

  void SliceMe(SubGrid *);       ///< Slice this grid according to the subgrid (internal function)

  Grid() = default;

 private:
  // Check if number is a power of 2
  bool isPow2(int);
  void makeDomainDecomposition();
};

/**
 * @brief SubGrid allows one to define a grid as a slice or averaged from another grid
 *
 */
class SubGrid {
 public:
  SubGrid(Grid * grid, SliceType type, int d, real x0):
          parentGrid(grid), type(type), direction(d) {
            idfx::pushRegion("SubGrid::SubGrid()");
            // Find the index of the current subgrid.
            auto x = Kokkos::create_mirror_view(parentGrid->x[direction]);
            Kokkos::deep_copy(x,parentGrid->x[direction]);
            int iref = -1;
            for(int i = 0 ; i < x.extent(0) - 1 ; i++) {
              if(x(i+1) >= x0) {
                if(x(i+1) - x0 > x0 - x(i) ) {
                  // we're closer to x(i)
                  iref = i;
                } else {
                  iref = i+1;
                }
                break;
              }
            }
            if(iref < 0) {
              std::stringstream msg;
              msg << "Cannot generate a subgrid in direction " << d << "with x0=" << x0
                  << "." << std::endl
                  << "Bounds are " << x(0) << "..." << x(x.extent(0)-1) << std::endl;

              IDEFIX_ERROR(msg);
            }
            this->index = iref;
            this->x0 = x(iref);

            this->grid = std::make_unique<Grid>(this);
            idfx::popRegion();
  }

  Grid *parentGrid;
  std::unique_ptr<Grid> grid;

  const SliceType type;
  const int direction;
  real x0;   ///< Cell center coordinate of the slice/average
  int index; ///< index in parent grid of the slice/average
};

#endif // GRID_HPP_
