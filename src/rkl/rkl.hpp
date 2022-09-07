// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef RKL_RKL_HPP_
#define RKL_RKL_HPP_

#include <vector>
#include "idefix.hpp"
#include "input.hpp"
#include "dataBlock.hpp"
#ifdef WITH_MPI
#include "mpi.hpp"
#endif


class RKLegendre {
 public:
  void Init(Input &, DataBlock &);
  void Cycle();
  void ResetStage();
  void ResetFlux();
  void EvolveStage(real);
  template <int> void CalcParabolicRHS(real);
  void ComputeDt();
  void ShowConfig();
  void Copy(IdefixArray4D<real>&, IdefixArray4D<real>&);

  IdefixArray4D<real> dU;      // variation of main cell-centered conservative variables
  IdefixArray4D<real> dU0;      // dU of the first stage
  IdefixArray4D<real> Uc0;      // Uc at initial stage
  IdefixArray4D<real> Uc1;      // Uc of the previous stage, Uc1 = Uc(stage-1)

  IdefixArray4D<real> dB;      // Variation of cell-centered magnetic variables
  IdefixArray4D<real> dB0;     // dB of the first stage
  IdefixArray4D<real> Vs0;     // Vs of initial stage
  IdefixArray4D<real> Vs1;     // Vs of previous stage

  #ifdef EVOLVE_VECTOR_POTENTIAL
  IdefixArray4D<real> dA;      // Variation of edge-centered vector potential
  IdefixArray4D<real> dA0;     // dA of the first stage
  IdefixArray4D<real> Ve0;     // Ve of initial stage
  IdefixArray4D<real> Ve1;     // Ve of previous stage
  #endif

  IdefixArray1D<int> varList;  // List of variables which should be evolved
  int nvarRKL{0};               // # of active variables

  real dt, cfl_rkl, rmax_par;
  int stage{0};

 private:
  void SetBoundaries(real);        // Enforce boundary conditions on the variables solved by RKL

  DataBlock *data;

#ifdef WITH_MPI
  Mpi mpi;                      // RKL-specific MPI layer
#endif

  bool haveVs{false};           // Whether we have (and need to compute) face-centered variables
  bool haveVc{false};           // Whether we need to compute cell-centered variables
  void AddVariable(int, std::vector<int> & );

  bool checkNan{false};         // whether we should look for Nans when RKL is running

 private:
  template<int> void LoopDir(real);   // Dimensional loop
};

#endif // RKL_RKL_HPP_
