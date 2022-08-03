// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_HDSOLVERS_HLLHD_HPP_
#define HYDRO_HDSOLVERS_HLLHD_HPP_

#include "../idefix.hpp"
#include "hydro.hpp"
#include "slopeLimiter.hpp"
#include "fluxHD.hpp"
#include "convertConsToPrimHD.hpp"

// Compute Riemann fluxes from states using HLL solver
template<const int DIR>
void Hydro::HllHD() {
  idfx::pushRegion("Hydro::HLL_Solver");

  int ioffset,joffset,koffset;
  ioffset=joffset=koffset=0;
  // Determine the offset along which we do the extrapolation
  if(DIR==IDIR) ioffset=1;
  if(DIR==JDIR) joffset=1;
  if(DIR==KDIR) koffset=1;

  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Vs = this->Vs;
  IdefixArray4D<real> Flux = this->FluxRiemann;
  IdefixArray3D<real> cMax = this->cMax;
  IdefixArray3D<real> csIsoArr = this->isoSoundSpeedArray;

  // Required for high order interpolations
  IdefixArray1D<real> dx = this->data->dx[DIR];

  [[maybe_unused]] real gamma = this->gamma;
  [[maybe_unused]] real gamma_m1 = gamma - ONE_F;
  [[maybe_unused]] real csIso = this->isoSoundSpeed;
  [[maybe_unused]] HydroModuleStatus haveIsoCs = this->haveIsoSoundSpeed;

  SlopeLimiter<DIR,NVAR> slopeLim(Vc,data->dx[DIR],shockFlattening);

  idefix_for("HLL_Kernel",
             data->beg[KDIR],data->end[KDIR]+koffset,
             data->beg[JDIR],data->end[JDIR]+joffset,
             data->beg[IDIR],data->end[IDIR]+ioffset,
    KOKKOS_LAMBDA (int k, int j, int i) {
      // Init the directions (should be in the kernel for proper optimisation by the compilers)
      constexpr int Xn = DIR+MX1;

      // Primitive variables
      real vL[NVAR];
      real vR[NVAR];

      // Conservative variables
      real uL[NVAR];
      real uR[NVAR];

      // Flux (left and right)
      real fluxL[NVAR];
      real fluxR[NVAR];

      // Signal speeds
      real cL, cR, cmax;

      // 1-- Store the primitive variables on the left, right, and averaged states
      slopeLim.ExtrapolatePrimVar(i, j, k, vL, vR);

      // 2-- Get the wave speed
#if HAVE_ENERGY
      cL = std::sqrt( gamma *(vL[PRS]/vL[RHO]) );
      cR = std::sqrt( gamma *(vR[PRS]/vR[RHO]) );
#else
      if(haveIsoCs == UserDefFunction) {
        cL = HALF_F*(csIsoArr(k,j,i)+csIsoArr(k-koffset,j-joffset,i-ioffset));
      } else {
        cL = csIso;
      }
      cR = cL;
#endif

      // 4.1
      real cminL = vL[Xn] - cL;
      real cmaxL = vL[Xn] + cL;

      real cminR = vR[Xn] - cR;
      real cmaxR = vR[Xn] + cR;

      real SL = FMIN(cminL, cminR);
      real SR = FMAX(cmaxL, cmaxR);

      cmax  = FMAX(FABS(SL), FABS(SR));

      // 2-- Compute the conservative variables: do this by extrapolation
      K_PrimToCons(uL, vL, gamma_m1);
      K_PrimToCons(uR, vR, gamma_m1);

      // 3-- Compute the left and right fluxes
      K_Flux(fluxL, vL, uL, cL*cL, Xn);
      K_Flux(fluxR, vR, uR, cR*cR, Xn);

      // 5-- Compute the flux from the left and right states
      if (SL > 0) {
#pragma unroll
        for (int nv = 0 ; nv < NFLX; nv++) {
          Flux(nv,k,j,i) = fluxL[nv];
        }
      } else if (SR < 0) {
#pragma unroll
        for (int nv = 0 ; nv < NFLX; nv++) {
          Flux(nv,k,j,i) = fluxR[nv];
        }
      } else {
#pragma unroll
        for(int nv = 0 ; nv < NFLX; nv++) {
          Flux(nv,k,j,i) = SL*SR*uR[nv] - SL*SR*uL[nv] + SR*fluxL[nv] - SL*fluxR[nv];
          Flux(nv,k,j,i) /= (SR - SL);
        }
      }

      //6-- Compute maximum wave speed for this sweep
      cMax(k,j,i) = cmax;
    }
  );

  idfx::popRegion();
}

#endif // HYDRO_HDSOLVERS_HLLHD_HPP_
