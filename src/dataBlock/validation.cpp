// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/dataBlock/validation.cpp
//
// Last modified : 10/2023
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2020 - 2023)
// - Soufiane Baghdadi <soufiane.baghdadi@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2020)
// - Clément Robert <clement.robert@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2021 - 2023)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "dataBlock.hpp"
#include "fluid.hpp"

int DataBlock::CheckNan() {
  idfx::pushRegion("DataBlock::Check");
  int nNans = hydro->CheckNan();
  if(haveDust) {
    for(int n = 0 ; n < dust.size() ; n++) {
      nNans += dust[n]->CheckNan();
    }
  }
  idfx::popRegion();
  return(nNans);
}

void DataBlock::Validate() {
  idfx::pushRegion("DataBlock::Validate");

  if(this->CheckNan()) {
    IDEFIX_ERROR("Nans were found in your initial conditions.");
  }

  idfx::popRegion();
}
