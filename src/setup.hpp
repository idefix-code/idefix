// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef SETUP_HPP_
#define SETUP_HPP_

#include "idefix.hpp"

class Setup {
 public:
  Setup();
  Setup(Input &, Grid &, DataBlock &);

  void InitFlow(DataBlock &);
  void MakeAnalysis(DataBlock&);
};

#endif // SETUP_HPP_
