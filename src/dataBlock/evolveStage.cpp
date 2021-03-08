// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "../idefix.hpp"
#include "dataBlock.hpp"

// Evolve one step forward in time of hydro
void DataBlock::EvolveStage() {
  idfx::pushRegion("DataBlock::EvolveStage");
  // Compute current when needed
  if(hydro.needCurrent) hydro.CalcCurrent();

  // Loop on all of the directions
  for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
    // Step one: extrapolate the variables to the sides, result is stored in the physics object
    hydro.ExtrapolatePrimVar(dir);

    // Step 2: compute the intercell flux with our Riemann solver, store the resulting InvDt
    hydro.CalcRiemannFlux(dir, this->t);

    // Step 2.5: compute intercell parabolic flux when needed
    if(hydro.haveParabolicTerms) hydro.CalcParabolicFlux(dir, this->t);

    // Step 3: compute the resulting evolution of the conserved variables, stored in Uc
    hydro.CalcRightHandSide(dir, this->t, this->dt);
  }

  // Step 4: add source terms to the conserved variables (curvature, rotation, etc)
  if(hydro.haveSourceTerms) hydro.AddSourceTerms(this->t, this->dt);

#if MHD == YES && DIMENSIONS >= 2
  // Compute the field evolution according to CT
  hydro.emf.CalcCornerEMF(this->t);
  if(hydro.haveResistivity || hydro.haveAmbipolar) hydro.emf.CalcNonidealEMF(this->t);
  hydro.emf.EnforceEMFBoundary();
  hydro.emf.EvolveMagField(this->t, this->dt);
  hydro.ReconstructVcField(hydro.Uc);
#endif

  idfx::popRegion();
}
