// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef GRAVITY_SELFGRAVITYITERATIVE_HPP_
#define GRAVITY_SELFGRAVITYITERATIVE_HPP_

#include <memory>
#include <vector>
#include <array>
#include "idefix.hpp"
#include "selfGravity.hpp"
#include "iterativesolver.hpp"
#include "laplacian.hpp"

template<typename Operator> class IterativeSolver;

class SelfGravityIterative final : public SelfGravity {
 public:
  void Init(Input &, DataBlock *) override;
  void ShowConfig() override;
  void InitSolver();
  void SolvePoisson() override;
  void AddSelfGravityPotential(IdefixArray3D<real> &) override;
  void EnrollUserDefBoundary(Laplacian::UserDefBoundaryFunc myFunc) override;
  void SubstractMeanDensity();

  std::unique_ptr<Laplacian> laplacian;
  std::unique_ptr<IterativeSolver<Laplacian>> iterativeSolver;


 private:
  real dt;  // CFL timestep

  std::array<Laplacian::LaplacianBoundaryType,3> lbound;
  std::array<Laplacian::LaplacianBoundaryType,3> rbound;
  bool havePreconditioner{false};
  GravitySolver solver{BICGSTAB};

  std::array<int,3> np_tot{0,0,0};
  IdefixArray3D<real> potential;
  IdefixArray3D<real> density;
  bool isPeriodic{true};
};

#endif // GRAVITY_SELFGRAVITYITERATIVE_HPP_
