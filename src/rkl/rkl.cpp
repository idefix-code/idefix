// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <cmath>
#include "rkl.hpp"
#include "dataBlock.hpp"


RKLegendre::RKLegendre() {
  // do nothing!
}


void RKLegendre::Init(DataBlock *datain) {
  idfx::pushRegion("RKLegendre::Init");

  // Save the datablock to which we are attached from now on
  this->data = datain;

  idfx::popRegion();
}


void RKLegendre::Cycle() {
  idfx::pushRegion("RKLegendre::Cycle");

  idfx::popRegion();
}

void RKLegendre::EvolveStage() {
  idfx::pushRegion("RKLegendre::EvolveStage");

  idfx::popRegion();
}


void RKLegendre::CalcParabolicRHS() {
  idfx::pushRegion("RKLegendre::CalcParabolicRHS");

  idfx::popRegion();
}
