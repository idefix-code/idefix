// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef DATABLOCK_BALANCEDSCHEME_HPP_
#define DATABLOCK_BALANCEDSCHEME_HPP_
#include "idefix.hpp"

class DataBlock;

class BalancedScheme {
 public:
  explicit BalancedScheme(DataBlock &data);
  void ComputeResidual(DataBlock &data);
  void LoadResidual(DataBlock &data);

 private:
  IdefixArray4D<real> dUc;
};

#endif // DATABLOCK_BALANCEDSCHEME_HPP_
