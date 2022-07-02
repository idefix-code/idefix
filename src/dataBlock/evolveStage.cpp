// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "../idefix.hpp"
#include "dataBlock.hpp"
#include "calcRightHandSide.hpp"
#include "calcParabolicFlux.hpp"
#include "calcRiemannFlux.hpp"

template<int dir> void DataBlock::LoopDir() {
    // Step 2: compute the intercell flux with our Riemann solver, store the resulting InvDt
    hydro.CalcRiemannFlux<dir>(this->t);

    // Step 2.5: compute intercell parabolic flux when needed
    if(hydro.haveExplicitParabolicTerms) hydro.CalcParabolicFlux<dir>(this->t);

    // Step 3: compute the resulting evolution of the conserved variables, stored in Uc
    hydro.CalcRightHandSide<dir>(this->t, this->dt);

    // Recursive: do next dimension
      LoopDir<dir+1>();
}

template<> void DataBlock::LoopDir<DIMENSIONS>() {
  // Do nothing
}


// Evolve one step forward in time of hydro
void DataBlock::EvolveStage() {
  idfx::pushRegion("DataBlock::EvolveStage");
  // Compute current when needed
  if(hydro.needExplicitCurrent) hydro.CalcCurrent();

  // enable shock flattening
  if(hydro.haveShockFlattening) hydro.shockFlattening.FindShock();

  // Loop on all of the directions
  LoopDir<IDIR>();

  // Step 4: add source terms to the conserved variables (curvature, rotation, etc)
  if(hydro.haveSourceTerms) hydro.AddSourceTerms(this->t, this->dt);

#if MHD == YES && DIMENSIONS >= 2
  // Compute the field evolution according to CT
  hydro.emf.CalcCornerEMF(this->t);
  if(hydro.resistivityStatus.isExplicit || hydro.ambipolarStatus.isExplicit) {
    hydro.emf.CalcNonidealEMF(this->t);
  }
  hydro.emf.EnforceEMFBoundary();
  #ifdef EVOLVE_VECTOR_POTENTIAL
    hydro.emf.EvolveVectorPotential(this->dt, hydro.Ve);
    hydro.emf.ComputeMagFieldFromA(hydro.Ve, hydro.Vs);
  #else
    hydro.emf.EvolveMagField(this->t, this->dt, hydro.Vs);
  #endif

  hydro.boundary.ReconstructVcField(hydro.Uc);
#endif

  idfx::popRegion();
}
