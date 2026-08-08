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
#include "laplacian.hpp"

#ifdef WITH_MPI
#include "mpi.hpp"
#endif

class DataBlock;

class SelfGravity {
 public:
  enum GravitySolver {JACOBI, BICGSTAB, PBICGSTAB, PCG, CG, PMINRES, MINRES, FFTSolver};

  virtual ~SelfGravity() = default;

  virtual void Init(Input &, DataBlock *) = 0;
  virtual void ShowConfig() = 0;
  virtual void SolvePoisson() = 0;
  virtual void AddSelfGravityPotential(IdefixArray3D<real> &) = 0;
  virtual void EnrollUserDefBoundary(Laplacian::UserDefBoundaryFunc myFunc) = 0;

  static std::unique_ptr<SelfGravity> Create(Input &, DataBlock *);

  real currentError{0};
  int nsteps{0};
  double elapsedTime{0.0};
  int skipSelfGravity{1};

 protected:
  DataBlock *data{nullptr};

};

#endif // GRAVITY_SELFGRAVITY_HPP_
