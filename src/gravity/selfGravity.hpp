// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef GRAVITY_SELFGRAVITY_HPP_
#define GRAVITY_SELFGRAVITY_HPP_

#include <vector>

#include "idefix.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "fluid_defs.hpp"
#include "iterativesolver.hpp"

#ifdef WITH_MPI
#include "mpi.hpp"
#endif

// Forward class hydro declaration
class DataBlock;

class SelfGravity {
 public:
  // Types of boundary which can be treated
  enum GravityBoundaryType {internalgrav, periodic, nullgrad, nullpot, userdef, axis, origin};
  enum GravitySolver {JACOBI, BICGSTAB, PBICGSTAB, PCG, CG, PMINRES, MINRES};

  // Handling userdef boundary.
  using UserDefBoundaryFunc = void (*) (DataBlock &, int dir, BoundarySide side,
                                       const real t, IdefixArray3D<real> &arr);

  SelfGravity();  // Default (empty) constructor
  void Init(Input &, DataBlock *);  // Initialisation of the class attributes
  void ShowConfig();                // display current configuration
  void InitSolver(); // (Re)initialisation of the solver for a given density distribution

  void InitInternalGrid(); // initialise the extra internal grid (for origin BCs)
  real ComputeJacobiCFL();  // Compute step for Jacobi following CFL conditions
  void SubstractMeanDensity();  // Compute and substract the average input density

  void ComputeLaplacian(IdefixArray3D<real> array,
                        IdefixArray3D<real> laplacian);
  // Compute the laplacian of the given array

  void EnforceBoundary(int dir, BoundarySide side, GravityBoundaryType type,
                       IdefixArray3D<real> &);
  // Enforce boundary conditions in a specific dir, side, following boundary type
  // and for the specified array

  void SetBoundaries(IdefixArray3D<real> &);  // Set the proper boundaries for the given array


  void SolvePoisson(); // Solve Poisson equation
  void AddSelfGravityPotential(IdefixArray3D<real> &);
  void EnrollUserDefBoundary(UserDefBoundaryFunc);  // Enroll user-defined boundary conditions
  void InitPreconditionner();   // For preconditionning versions

  // real gravCst; // 4*pi*G

  // User defined Boundary conditions
  UserDefBoundaryFunc userDefBoundaryFunc{NULL};
  bool haveUserDefBoundary{false};

  IterativeSolver<SelfGravity> *iterativeSolver;

  #ifdef WITH_MPI
  Mpi mpi;  // Mpi object when WITH_MPI is set
  #endif

  #ifdef DEBUG_GRAVITY
  // Used to get fields usefull for debugging
  void WriteField(std::ofstream &, IdefixArray3D<real> &, int = 0);
  void WriteField(std::ofstream &, IdefixArray1D<real> &, int = 0);
  std::ofstream rhoFile;
  std::ofstream potentialFile;
  std::ofstream geometryFile;
  #endif

  real currentError{0};       // last error of the iterative solver
  int nsteps{0};              // # of steps of the latest iteration
  double elapsedTime;        // time spent solving self gravity

  // Whether we should skip self-gravity computation every n steps
  int skipSelfGravity{1};

 private:
  DataBlock *data;  // My parent data object
  IdefixArray3D<real> potential;  // Gravitational potential
  IdefixArray3D<real> density;  // Density
  real dt;  // CFL timestep

  // Local potential array size
  std::vector<int> np_tot;
  std::vector<int> np_int;

  std::vector<int> nghost;
  std::vector<int> beg;
  std::vector<int> end;

  // offset in the left and right direction between selfgravity grid and datablock grid
  std::vector<int> loffset;
  std::vector<int> roffset;

  // Grid for self-gravity solver (note that this grid may extend the grid of the current datablock)
  std::vector<IdefixArray1D<real>> x;    ///< geometrical central points
  std::vector<IdefixArray1D<real>> dx;   ///< cell width

  IdefixArray1D<real> sinx2;            ///< sinx2 (only in spherical)

  IdefixArray3D<real> dV;                ///< cell volume
  std::vector<IdefixArray3D<real>> A;    ///< cell left interface area

  #ifdef WITH_MPI
  MPI_Comm originComm;                  ///< MPI communicator used by the origin boundary condition
  #endif

  bool isPeriodic; // Periodicity status of the density distribution
  GravityBoundaryType lbound[3];  // Boundary condition to the left
  GravityBoundaryType rbound[3];  // Boundary condition to the right
                           // Warning : might differ from (M)HD solver !

  IdefixArray3D<real> precond; // Diagonal preconditionner

  #ifdef WITH_MPI
  IdefixArray4D<real> arr4D; // Intermediate array for boundary handling
  #endif

  bool isTwoPi{false};

  GravitySolver solver; // The solver  used to solve Poisson

  bool havePreconditioner{false}; // Use of preconditionner (or not)
};

#endif // GRAVITY_SELFGRAVITY_HPP_
