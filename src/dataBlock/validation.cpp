// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "dataBlock.hpp"
#include "fluid.hpp"

int DataBlock::CheckNan() {
  idfx::pushRegion("DataBlock::Validate");
  int nNans = hydro->CheckNan();
  if(haveDust) {
    for(int n = 0 ; n < dust.size() ; n++) {
      nNans += dust[n]->CheckNan();
    }
  }
  return(nNans);
}

void DataBlock::Validate() {
  idfx::pushRegion("DataBlock::Validate");

  if(this->CheckNan()) {
    IDEFIX_ERROR("Nans were found in your initial conditions.");
  }

  idfx::popRegion();
}
