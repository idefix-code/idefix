// ***********************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************************

#include "../idefix.hpp"
#include "hydro.hpp"

// Compute parabolic fluxes
void Hydro::CalcParabolicFlux(int dir, const real t) {
  idfx::pushRegion("Hydro::CalcParabolicFlux");

  if(this->haveResistivity || this->haveAmbipolar) {
    this->AddNonIdealMHDFlux(dir,t);
  }

  idfx::popRegion();  
}
