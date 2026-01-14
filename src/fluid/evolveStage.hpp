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
#include "forcing.hpp"

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
  if(haveDrag) {
    if(!drag->IsImplicit()) {
      drag->AddDragForce(dt);
    }
  }

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


// Evolve one step forward in time of forcing
template<typename Phys>
void Fluid<Phys>::EvolveForcing(const real t, const real dt) {
  idfx::pushRegion("Fluid::EvolveForcing");

  data->forcing->ComputeForcing(dt);
  data->forcing->oUprocesses.WriteTimestep(t, dt);
  if (data->forcing->write) {
    data->forcing->oUprocesses.WriteProcessesValues(t);
    data->forcing->oUprocesses.WriteNormalValues(t);
  }

  if (data->forcing->stillHaveForcing and t > data->forcing->stopTime) data->forcing->stillHaveForcing = false;

  // Loop on all of the directions
  LoopForcingDir<IDIR>(t,dt);

  idfx::popRegion();
}

template<typename Phys>
template<int dir>
void Fluid<Phys>::LoopForcingDir(const real t, const real dt) {

    CalcForcingRHS<dir>(t,dt);

    // Recursive: do next dimension
    if constexpr (dir+1 < DIMENSIONS) LoopForcingDir<dir+1>(t, dt);
}

template<typename Phys, int dir>
struct Fluid_CalcForcingFunctor {
  //*****************************************************************
  // Functor constructor
  //*****************************************************************
  explicit Fluid_CalcForcingFunctor (Fluid<Phys> *hydro, real dt) {
    Uc   = hydro->Uc;
    Vc   = hydro->Vc;
    Flux = hydro->FluxRiemann;
    dV   = hydro->data->dV;
    forcingTerm = hydro->data->forcing->forcingTerm;

    rt   = hydro->data->rt;
    dmu  = hydro->data->dmu;
    dx   = hydro->data->dx[dir];
    dx2  = hydro->data->dx[JDIR];
//    invDt = hydro->InvDt;
//    cMax = hydro->cMax;
//    dMax = hydro->dMax;
    this->dt = dt;
  }
  //*****************************************************************
  // Functor Variables
  //*****************************************************************
  IdefixArray4D<real> Uc;
  IdefixArray4D<real> Vc;
  IdefixArray3D<real> dV;
  IdefixArray4D<real> Flux;
  IdefixArray4D<real> forcingTerm;

  IdefixArray1D<real> rt;
  IdefixArray1D<real> dmu;
  IdefixArray1D<real> dx;
  IdefixArray1D<real> dx2;
  real dt;

  //*****************************************************************
  // Functor Operator
  //*****************************************************************
  KOKKOS_INLINE_FUNCTION void operator() (const int k, const int j,  const int i) const {
    const int ioffset = (dir==IDIR) ? 1 : 0;
    const int joffset = (dir==JDIR) ? 1 : 0;
    const int koffset = (dir==KDIR) ? 1 : 0;

    real dtdV=dt / dV(k,j,i);
    real rhs[Phys::nvar];

    #pragma unroll
    for(int nv = 0 ; nv < Phys::nvar ; nv++) {
      rhs[nv] = ZERO_F;
    }


    // elmentary length for gradient computations
    const int ig = ioffset*i + joffset*j + koffset*k;
    real dl = dx(ig);
    #if GEOMETRY == POLAR
      if constexpr (dir==JDIR)
        dl = dl*x1(i);

    #elif GEOMETRY == SPHERICAL
      if constexpr(dir==JDIR)
        dl = dl*rt(i);
      else if constexpr(dir==KDIR)
          dl = dl*rt(i)*dmu(j)/dx2(j);
    #endif


      rhs[MX1+dir] = dt * Vc(RHO,k,j,i) * forcingTerm(dir,k,j,i);
      if constexpr(Phys::pressure) {
        //  rho * v . f
        rhs[ENG] += dt * Vc(RHO,k,j,i) * Vc(VX1+dir,k,j,i) * forcingTerm(dir,k,j,i);
      } // Pressure

      // Particular cases if we do not sweep all of the components
      #if DIMENSIONS == 1 && COMPONENTS > 1
        EXPAND(                                                           ,
                  rhs[MX2] += dt * Vc(RHO,k,j,i) * forcingTerm(JDIR,k,j,i);   ,
                  rhs[MX3] += dt * Vc(RHO,k,j,i) * forcingTerm(KDIR,k,j,i);    )
        if constexpr(Phys::pressure) {
          rhs[ENG] += dt * (EXPAND( ZERO_F                                              ,
                                    + Vc(RHO,k,j,i) * Vc(VX2,k,j,i) * forcingTerm(JDIR,k,j,i)   ,
                                    + Vc(HO,k,j,i) * Vc(VX3,k,j,i) * forcingTerm(KDIR,k,j,i) ));
        }
      #endif
      #if DIMENSIONS == 2 && COMPONENTS == 3
        // Only add this term once!
        if constexpr (dir==JDIR) {
          rhs[MX3] += dt * Vc(RHO,k,j,i) * forcingTerm(KDIR,k,j,i);
          if constexpr(Phys::pressure) {
            rhs[ENG] += dt * Vc(RHO,k,j,i) * Vc(VX3,k,j,i) * forcingTerm(KDIR,k,j,i);
          }
        }
      #endif

  #pragma unroll
  for(int nv = 0 ; nv < Phys::nvar ; nv++) {
    Uc(nv,k,j,i) = Uc(nv,k,j,i) + rhs[nv];
  }

  // WARNING DO SOMETHING FOR THE TIMESTEP DT?

  }
};

template<typename Phys>
template<int dir>
void Fluid<Phys>::CalcForcingRHS(const real t, const real dt) {
  idfx::pushRegion("Fluid::CalcForcing");

  auto calcForcing = Fluid_CalcForcingFunctor<Phys,dir>(this,dt);
  /////////////////////////////////////////////////////////////////////////////
  // Final conserved quantity budget from fluxes divergence
  /////////////////////////////////////////////////////////////////////////////
  idefix_for("CalcForcing",
             data->beg[KDIR],data->end[KDIR],
             data->beg[JDIR],data->end[JDIR],
             data->beg[IDIR],data->end[IDIR],
              calcForcing);

  idfx::popRegion();
}

#endif //FLUID_EVOLVESTAGE_HPP_
