// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/dataBlock/evolveStage.cpp
//
// Last modified : 02/2025
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2020 - 2025)
// - Soufiane Baghdadi <soufiane.baghdadi@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2020)
// - Clément Robert <clement.robert@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2021)
// and other code contributors
//
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
    // Add implicit term for dust drag
    if(dust[0]->drag->IsImplicit()) {
      for(int i = 0 ; i < dust.size() ; i++) {
        dust[i]->drag->AddImplicitBackReaction(this->dt,dust[0]->drag->implicitFactor);
      }
      dust[0]->drag->NormalizeImplicitBackReaction(this->dt);
      for(int i = 0 ; i < dust.size() ; i++) {
        dust[i]->drag->AddImplicitFluidMomentum(this->dt);
      }
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
