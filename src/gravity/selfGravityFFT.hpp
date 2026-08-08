// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef GRAVITY_SELFGRAVITYFFT_HPP_
#define GRAVITY_SELFGRAVITYFFT_HPP_

#include <memory>
#include <array>
#include <KokkosFFT.hpp>

#include "fft.hpp"
#include "selfGravity.hpp"


class SelfGravityFFT final : public SelfGravity {
 public:
  void Init(Input &, DataBlock *) override;
  void ShowConfig() override;
  void SolvePoisson() override;
  void AddSelfGravityPotential(IdefixArray3D<real> &) override;
  void EnrollUserDefBoundary(Laplacian::UserDefBoundaryFunc myFunc) override;
  void SetBoundaries(IdefixArray3D<real> &);
  void EnforcePeriodic(int dir, BoundarySide side, IdefixArray3D<real> &);
 private:
  void SubstractMeanDensity();
  void CheckCompatibility();

  std::unique_ptr<FFT> fft;  ///< FFT wrapper

  // FFT work arrays on active domain (no ghosts)
  IdefixArray3D<real> rho;
  IdefixArray3D<real> phi;
  IdefixArray3D<complex> rhoF;
  IdefixArray3D<complex> phiF;
  std::array<int,3> npf{1,1,1}; // [k,j,i]
  std::array<int,3> npf_t{1,1,1}; // [j,k,i]
  std::array<int,3> npr{1,1,1}; // [k,j,i]
  std::array<int,3> npr_glob{1,1,1}; // [k,j,i]

  std::array<IdefixArray1D<real>,3> kx_glob;
  std::array<IdefixArray1D<real>,3> kx;

  std::array<int,3> begin{0,0,0}; // [k,j,i]
  std::array<int,3> end{0,0,0}; // [k,j,i]
};

#endif // GRAVITY_SELFGRAVITYFFT_HPP_
