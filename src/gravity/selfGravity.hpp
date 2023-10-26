// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef GRAVITY_SELFGRAVITY_HPP_
#define GRAVITY_SELFGRAVITY_HPP_

#include <memory>
#include <vector>

#include "idefix.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "fluid_defs.hpp"
#include "iterativesolver.hpp"
#include "laplacian.hpp"

#ifdef WITH_MPI
#include "mpi.hpp"
#endif

// Forward class hydro declaration
class DataBlock;

class SelfGravity {
 public:
  enum GravitySolver {JACOBI, BICGSTAB, PBICGSTAB, PCG, CG, PMINRES, MINRES};

  void Init(Input &, DataBlock *);  // Initialisation of the class attributes
  void ShowConfig();                // display current configuration
  void InitSolver(); // (Re)initialisation of the solver for a given density distribution

  void SubstractMeanDensity();  // Compute and substract the average input density

  void SolvePoisson(); // Solve Poisson equation
  void AddSelfGravityPotential(IdefixArray3D<real> &);

  void EnrollUserDefBoundary(Laplacian::UserDefBoundaryFunc myFunc);  // User-defined boundary

  IterativeSolver<Laplacian> *iterativeSolver;

  // The linear operator involved in Poisson equation
  std::unique_ptr<Laplacian> laplacian;

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
  std::array<int,3> np_tot;

  std::array<Laplacian::LaplacianBoundaryType,3> lbound;  // Boundary condition to the left
  std::array<Laplacian::LaplacianBoundaryType,3> rbound;  // Boundary condition to the right

  bool isPeriodic;
  bool havePreconditioner{false};
  GravitySolver solver; // The solver  used to solve Poisson
};

#endif // GRAVITY_SELFGRAVITY_HPP_
