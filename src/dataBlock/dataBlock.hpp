// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef DATABLOCK_DATABLOCK_HPP_
#define DATABLOCK_DATABLOCK_HPP_

#include <vector>
#include <string>
#include <map>
#include <memory>

#include "idefix.hpp"
#include "grid.hpp"
#include "gridHost.hpp"
#include "output.hpp"
#include "gravity.hpp"
#include "stateContainer.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////////
/// The DataBlock class is designed to store the data and child class instances that belongs to the
/// current MPI process ONLY.  In particular grid-related arrays of the DataBlock
/// only contains information of the local portion of the grid that belongs to the current MPI
/// process. The full grid is defined in the parent Grid class (for which a pointer is defined).
/// Note that all of the arrays of a DataBlock are on the device. If a Host access is needed, it
/// is recommended to use a DataBlockHost instance from this DataBlock, and sync it.
/////////////////////////////////////////////////////////////////////////////////////////////////

// forward class declaration (used by enrollment functions)
class DataBlock;
class Vtk;
class Dump;
class Fargo;
class Gravity;
class PlanetarySystem;
template<typename Phys>
class Fluid;

// Forward class hydro declaration
#include "physics.hpp"

using GridCoarseningFunc = void(*) (DataBlock &);

using StepFunc = void (*) (DataBlock &, const real t, const real dt);

class DataBlock {
 public:
  // Local grid information
  std::vector<IdefixArray1D<real>> x;    ///< geometrical central points
  std::vector<IdefixArray1D<real>> xr;   ///< cell right interface
  std::vector<IdefixArray1D<real>> xl;   ///< cell left interface
  std::vector<IdefixArray1D<real>> dx;   ///< cell width
  std::vector<IdefixArray1D<real>> xgc;  ///< cell geometrical cell center
  IdefixArray1D<real> rt;      ///< In spherical coordinates, gives $\tilde{r}$
  IdefixArray1D<real> sinx2m;  ///< In spherical coordinates,
                               ///< gives sin(th) at a j-1/2 interface
  IdefixArray1D<real> tanx2m;  ///< In spherical coordinates,
                               ///< gives tan(th) at a j-1/2 interface
  IdefixArray1D<real> sinx2;   ///< In spherical coordinates, gives sin(th) at the cell center
  IdefixArray1D<real> tanx2;   ///< In spherical coordinates, gives tan(th) at the cell center
  IdefixArray1D<real> dmu;     ///< In spherical coordinates,
                               ///< gives the $\theta$ volume = fabs(cos(th_m) - cos(th_p))

  std::vector<IdefixArray2D<int>> coarseningLevel; ///< Grid coarsening levels
                                                  ///< (only defined when coarsening
                                                  ///< is enabled)
  std::vector<bool> coarseningDirection;  ///< whether a coarsening is used in each direction

  std::vector<real> xbeg;             ///< Beginning of active domain in datablock
  std::vector<real> xend;             ///< End of active domain in datablock

  IdefixArray3D<real> dV;                ///< cell volume
  std::vector<IdefixArray3D<real>> A;    ///< cell left interface area

  std::vector<int> np_tot;        ///< total number of grid points in datablock
  std::vector<int> np_int;        ///< active number of grid points in datablock (excl. ghost cells)

  std::vector<int> nghost;               ///< number of ghost cells at each boundary
  std::vector<BoundaryType> lbound;      ///< Boundary condition to the left
  std::vector<BoundaryType> rbound;      ///< Boundary condition to the right

  bool haveAxis{false};       ///< DataBlock contains points on the axis and a special treatment
                               ///< has been required for these.

  GridCoarsening haveGridCoarsening{GridCoarsening::disabled};
                                ///< Is grid coarsening enabled?
  GridCoarseningFunc gridCoarseningFunc{NULL};
                               ///< The user-defined grid coarsening level computation function



  std::vector<int> beg;       ///< First local index of the active domain
  std::vector<int> end;       ///< Last local index of the active domain+1

  std::vector<int> gbeg;      ///< First global index of the active domain of this datablock
  std::vector<int> gend;      ///< Last global index of the active domain of this datablock

  real dt;                     ///< Current timestep
  real t;                      ///< Current time

  Grid *mygrid;                ///< Parent grid object

  std::map<std::string, StateContainer> states;
                                ///< conservative state of the datablock
                                ///< (contains references to dedicated objects)

  std::unique_ptr<Fluid<DefaultPhysics>> hydro;   ///< The Hydro object attached to this datablock
  std::vector<std::unique_ptr<Fluid<DustPhysics>>> dust; ///< Holder for zero pressure dust fluid

  std::unique_ptr<Vtk> vtk;
  std::unique_ptr<Dump> dump;

  DataBlock(Grid &, Input &);     ///< init from a Grid object

  void MakeGeometry();                ///< Compute geometrical terms
  void DumpToFile(std::string);   ///< Dump current datablock to a file for inspection
  void Validate();                ///< error out early in case problems are found in IC
  int CheckNan();                 ///< Return the number of cells which have Nans

  // The Planetary system
  bool haveplanetarySystem{false};
  std::unique_ptr<PlanetarySystem> planetarySystem;


  bool rklCycle{false};           ///<  // Set to true when we're inside a RKL call

  void EvolveStage();             ///< Evolve this DataBlock by dt
  void EvolveRKLStage();          ///< Evolve this DataBlock by dt for terms impacted by RKL
  void SetBoundaries();       ///< Enforce boundary conditions to this datablock
  void DeriveVectorPotential(); ///< Compute magnetic fields from vector potential where applicable
  void Coarsen();             ///< Coarsen this datablock and its objects
  void ShowConfig();              ///< Show the datablock's configuration
  real ComputeTimestep();         ///< compute maximum timestep from current state of affairs

  void ResetStage();              ///< Reset the variables needed at each major integration Stage

  void EnrollGridCoarseningLevels(GridCoarseningFunc);
                                  ///< Enroll a user function to compute coarsening levels
  void CheckCoarseningLevels();   ///< Check that coarsening levels satisfy requirements

  // Do we use fargo-like scheme ? (orbital advection)
  bool haveFargo{false};
  std::unique_ptr<Fargo> fargo;

  // Do we have Gravity ?
  bool haveGravity{false};
  std::unique_ptr<Gravity> gravity;

  // User step functions (before or after the main integrator step)
  void LaunchUserStepFirst();     ///< perform user-defined step before main integration step
  void LaunchUserStepLast();      ///< Perform user-defined step after main integration step

  void EnrollUserStepFirst(StepFunc);
  void EnrollUserStepLast(StepFunc);

 private:
  void WriteVariable(FILE* , int , int *, char *, void*);
  void ComputeGridCoarseningLevels();   ///< Call user defined function to define Coarsening levels

  // User Steps (either before or after the main integration loop)
  bool haveUserStepFirst{false};
  bool haveUserStepLast{false};

  StepFunc userStepFirst{nullptr};
  StepFunc userStepLast{nullptr};
};

#endif // DATABLOCK_DATABLOCK_HPP_
