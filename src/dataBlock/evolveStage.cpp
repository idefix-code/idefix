// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "../idefix.hpp"
#include "dataBlock.hpp"
#include "fluid.hpp"

// Evolve one step forward in time of hydro
void DataBlock::EvolveStage() {
  idfx::pushRegion("DataBlock::EvolveStage");

  hydro->EvolveStage(this->t,this->dt);

  if(haveDust) {
    for(int i = 0 ; i < dust.size() ; i++) {
      dust[i]->EvolveStage(this->t,this->dt);
    }
  }

  idfx::popRegion();
}

void DataBlock::EvolveRKLStage() {
  idfx::pushRegion("DataBlock::EvolveRKLStage");
  if(hydro->haveRKLParabolicTerms) {
    hydro->rkl->Cycle();
  }
  idfx::popRegion();
}
