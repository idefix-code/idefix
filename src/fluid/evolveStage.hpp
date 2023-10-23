// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_EVOLVESTAGE_HPP_
#define FLUID_EVOLVESTAGE_HPP_

#include "fluid.hpp"
#include "riemannSolver.hpp"
template<typename Phys>
template<int dir>
void Fluid<Phys>::LoopDir(const real t, const real dt) {
    // Step 2: compute the intercell flux with our Riemann solver, store the resulting InvDt
    this->rSolver->template CalcFlux<dir>(this->FluxRiemann);

    // Step 2.5: compute intercell parabolic flux when needed
    if(haveExplicitParabolicTerms) CalcParabolicFlux<dir>(t);

    // If we have tracers, compute the tracer intercell flux
    if(haveTracer) {
      this->tracer->template CalcFlux<dir, Phys>(this->FluxRiemann);
    }

    // Step 3: compute the resulting evolution of the conserved variables, stored in Uc
    CalcRightHandSide<dir>(t,dt);
    if(haveTracer) {
      this->tracer->template CalcRightHandSide<dir, Phys>(this->FluxRiemann,t ,dt);
    }

    // Recursive: do next dimension
    if constexpr (dir+1 < DIMENSIONS) LoopDir<dir+1>(t, dt);
}



// Evolve one step forward in time of hydro
template<typename Phys>
void Fluid<Phys>::EvolveStage(const real t, const real dt) {
  idfx::pushRegion("Fluid::EvolveStage");
  // Compute current when needed
  if(needExplicitCurrent) CalcCurrent();

  if(hallStatus.status == UserDefFunction) {
    if(hallDiffusivityFunc)
      hallDiffusivityFunc(*data, t, xHall);
    else
      IDEFIX_ERROR("No user-defined Hall diffusivity function has been enrolled");
  }

  if constexpr(Phys::eos) {
    eos->Refresh(*data, t);
  }

  // Loop on all of the directions
  LoopDir<IDIR>(t,dt);

  // Step 4: add source terms to the conserved variables (curvature, rotation, etc)
  if(haveSourceTerms) AddSourceTerms(t, dt);

  // Step 5: add drag when needed
  if(haveDrag) drag->AddDragForce(dt);

  if constexpr(Phys::mhd) {
    #if DIMENSIONS >= 2
      // Compute the field evolution according to CT
      emf->CalcCornerEMF(t);
      if(resistivityStatus.isExplicit || ambipolarStatus.isExplicit) {
        emf->CalcNonidealEMF(t);
      }
      emf->EnforceEMFBoundary();
      #ifdef EVOLVE_VECTOR_POTENTIAL
        emf->EvolveVectorPotential(dt, Ve);
        emf->ComputeMagFieldFromA(Ve, Vs);
      #else
        emf->EvolveMagField(t, dt, Vs);
      #endif

      boundary->ReconstructVcField(Uc);
    #endif
  }

  idfx::popRegion();
}
#endif //FLUID_EVOLVESTAGE_HPP_
