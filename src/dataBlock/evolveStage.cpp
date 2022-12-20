// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "../idefix.hpp"
#include "dataBlock.hpp"
#include "calcRightHandSide.hpp"
#include "calcParabolicFlux.hpp"
#include "calcRiemannFlux.hpp"

// Evolve one step forward in time of hydro
void DataBlock::EvolveStage() {
  idfx::pushRegion("DataBlock::EvolveStage");
  
  hydro->EvolveStage(this->t,this->dt);

  idfx::popRegion();
}
