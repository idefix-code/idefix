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
  void Copy(IdefixArray4D<real>&, IdefixArray4D<real>&);

  IdefixArray4D<real> dU;      // variation of main cell-centered conservative variables
  IdefixArray4D<real> dU0;      // dU of the first stage
  IdefixArray4D<real> Uc1;      // Uc of the previous stage, Uc1 = Uc(stage-1)

  IdefixArray4D<real> dB;      // Variation of cell-centered magnetic variables
  IdefixArray4D<real> dB0;     // dB of the first stage
  IdefixArray4D<real> Vs1;     // Vs of previous stage

  IdefixArray1D<int> varList;  // List of variables which should be evolved
  int nvarRKL{0};               // # of active variables

  real dt, cfl_rkl, rmax_par;
  int stage{0};

 private:
  DataBlock *data;

  bool haveVs{false};           // Whether we have (and need to compute) cell-centered variables
  void AddVariable(int, IdefixArray1D<int>::HostMirror & );
};

#endif // RKL_RKL_HPP_
