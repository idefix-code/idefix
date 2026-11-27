// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_RIEMANNSOLVER_RADIATIONSOLVERS_LFRRAD_HPP_
#define FLUID_RIEMANNSOLVER_RADIATIONSOLVERS_LFRRAD_HPP_

#include "../idefix.hpp"
#include "fluid.hpp"
#include "extrapolateToFaces.hpp"
#include "flux.hpp"
#include "convertConsToPrim.hpp"
#include "speedRad.hpp"
#include "lim_fluxRad.hpp"

// Compute Riemann fluxes from states using Lax-Friedrichs-Rusanov solver
template <typename Phys>
template<const int DIR>
void RiemannSolver<Phys>::LFRRad(IdefixArray4D<real> &Flux) {
  idfx::pushRegion("RiemannSolver::LFR_Rad");

  constexpr int ioffset = (DIR==IDIR) ? 1 : 0;
  constexpr int joffset = (DIR==JDIR) ? 1 : 0;
  constexpr int koffset = (DIR==KDIR) ? 1 : 0;

  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray3D<real> cMax = this->cMax;

  // Required for high order interpolations
  IdefixArray1D<real> dx = this->data->dx[DIR];


  ExtrapolateToFaces<Phys,DIR> extrapol = *this->GetExtrapolator<DIR>();

  // Reduced velocity of light
  real reduced_c = this->reduced_c;

  // Type of rad flux limiter
  Type_Radlimiter rad_limiter = radLimiter;

  idefix_for("LFR_Rad_Kernel",
             data->beg[KDIR],data->end[KDIR]+koffset,
             data->beg[JDIR],data->end[JDIR]+joffset,
             data->beg[IDIR],data->end[IDIR]+ioffset,
    KOKKOS_LAMBDA (int k, int j, int i) {
      // Init the directions (should be in the kernel for proper optimisation by the compilers)
      constexpr int Xn = DIR+MX1;

      // Primitive variables
      real vL[Phys::nvar];
      real vR[Phys::nvar];
      real v[Phys::nvar];
      real voffset[Phys::nvar];

      for(int nv = 0 ; nv < Phys::nvar; nv++) {
        v[nv] = Vc(nv,k,j,i);
        voffset[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset);
      }

      // Conservative variables
      real uL[Phys::nvar];
      real uR[Phys::nvar];

      // Flux (left and right)
      real fluxL[Phys::nvar];
      real fluxR[Phys::nvar];

      //VWave speeds
      real lambdaL[2];
      real lambdaR[2];

      // xi from closure
      real xiL, xiR;

      // 1-- Store the primitive variables on the left, right, and averaged states
      extrapol.ExtrapolatePrimVar(i, j, k, vL, vR);

      // Limit the fluxes after extrapolation to satisfy Fr<=Er
      K_LimitRadFlux<Phys>(vL,vR,v,voffset,rad_limiter);

      // 2-- Get the wave speed
      K_SpeedsRad(lambdaL,vL,Xn, reduced_c,&xiL);
      K_SpeedsRad(lambdaR,vR,Xn, reduced_c,&xiR);

      real lambda_max_L = FMAX(lambdaL[0],lambdaL[1]);
      real lambda_max_R = FMAX(lambdaR[0],lambdaR[1]);
      real lambda_min_L = FMIN(lambdaL[0],lambdaL[1]);
      real lambda_min_R = FMIN(lambdaR[0],lambdaR[1]);

      real SR = FMAX(lambda_max_L,lambda_max_R);
      real SL = FMIN(lambda_min_L,lambda_min_R);

      real cmax  = FABS(FMAX(SL, SR));

      // 3-- Compute the conservative variables: do this by extrapolation
      K_PrimToCons<Phys>(uL, vL, NULL);
      K_PrimToCons<Phys>(uR, vR, NULL);

      // 4-- Compute the left and right fluxes (wave speed is null)
      K_Flux<Phys,DIR>(fluxL, vL, uL, reduced_c,xiL);
      K_Flux<Phys,DIR>(fluxR, vR, uR, reduced_c,xiR);

      // 5-- Compute the flux from the left and right states
#pragma unroll
      for(int nv = 0 ; nv < Phys::nvar; nv++) {
        Flux(nv,k,j,i) = 0.5*(fluxL[nv] + fluxR[nv]-cmax*(uR[nv]-uL[nv]));
      }

      //6-- Compute maximum wave speed for this sweep
      cMax(k,j,i) = cmax;
    }
  );

  idfx::popRegion();
}

#endif // FLUID_RIEMANNSOLVER_RADIATIONSOLVERS_LFRRAD_HPP_
