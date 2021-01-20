// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef SETUP_HPP_
#define SETUP_HPP_

#include "idefix.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "dataBlock.hpp"

// These two will likely be used in the code, so we include them here
#include "gridHost.hpp"
#include "dataBlockHost.hpp"

class Setup {
 public:
  Setup();
  Setup(Input &, Grid &, DataBlock &);

  void InitFlow(DataBlock &);
  void MakeAnalysis(DataBlock&);
};

#endif // SETUP_HPP_
