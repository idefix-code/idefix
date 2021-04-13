// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef RKL_RKL_HPP_
#define RKL_RKL_HPP_

#include "idefix.hpp"

// forward class declaration
class DataBlock;

class RKLegendre {
 public:
  RKLegendre();
  void Init(DataBlock *);
  void Cycle();
  void EvolveStage();
  void CalcParabolicRHS();

 private:
  DataBlock *data;
  int stage;
};

#endif // RKL_RKL_HPP_
