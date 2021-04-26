// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef RKL_RKL_HPP_
#define RKL_RKL_HPP_

#include "idefix.hpp"
#include "input.hpp"
#include "dataBlock.hpp"


class RKLegendre {
 public:
  RKLegendre();
  void Init(Input &, DataBlock &);
  void Cycle();
  void ResetStage();
  void ResetFlux();
  void EvolveStage(real);
  void CalcParabolicRHS(int, real);
  void ComputeDt();

  IdefixArray4D<real> dU;      // variation of main cell-centered conservative variables
  IdefixArray4D<real> dU0;      // dU of the first stage
  IdefixArray4D<real> Uc1;      // Uc of the previous stage, Uc1 = Uc(stage-1)

  real dt, cfl_rkl, rmax_par;

 private:
  DataBlock *data;
  int stage;
};

#endif // RKL_RKL_HPP_
