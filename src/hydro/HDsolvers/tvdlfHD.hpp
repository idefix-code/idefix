// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_HDSOLVERS_TVDLFHD_HPP_
#define HYDRO_HDSOLVERS_TVDLFHD_HPP_

#include "../idefix.hpp"

// Compute Riemann fluxes from states using TVDLF solver
template<const int DIR, const int Xn, const int Xt, const int Xb>
void Hydro::TvdlfHD() {
  idfx::pushRegion("Hydro::TVDLF_Solver");

  int ioffset,joffset,koffset;
  ioffset=joffset=koffset=0;
  // Determine the offset along which we do the extrapolation
  if(DIR==IDIR) ioffset=1;
  if(DIR==JDIR) joffset=1;
  if(DIR==KDIR) koffset=1;

  IdefixArray4D<real> PrimL = this->PrimL;
  IdefixArray4D<real> PrimR = this->PrimR;
  IdefixArray4D<real> Flux = this->FluxRiemann;
  IdefixArray3D<real> cMax = this->cMax;
  IdefixArray3D<real> csIsoArr = this->isoSoundSpeedArray;

  real gamma = this->gamma;
  real gamma_m1=this->gamma-ONE_F;
  real csIso = this->isoSoundSpeed;
  HydroModuleStatus haveIsoCs = this->haveIsoSoundSpeed;

  idefix_for("TVDLF_Kernel",
             data->beg[KDIR],data->end[KDIR]+koffset,
             data->beg[JDIR],data->end[JDIR]+joffset,
             data->beg[IDIR],data->end[IDIR]+ioffset,
    KOKKOS_LAMBDA (int k, int j, int i) {
      // Primitive variables
      real vL[NVAR];
      real vR[NVAR];
      real vRL[NVAR];

      // Conservative variables
      real uL[NVAR];
      real uR[NVAR];

      // Flux (left and right)
      real fluxL[NVAR];
      real fluxR[NVAR];

      // Signal speeds
      real cRL, cmax;

      // 1-- Read primitive variables
#pragma unroll
      for(int nv = 0 ; nv < NVAR; nv++) {
        vL[nv] = PrimL(nv,k,j,i);
        vR[nv] = PrimR(nv,k,j,i);
        vRL[nv] = HALF_F*(vL[nv]+vR[nv]);
      }

      // 2-- Get the wave speed
#if HAVE_ENERGY
      cRL = std::sqrt( (gamma_m1+ONE_F)*(vRL[PRS]/vRL[RHO]) );
#else
      if(haveIsoCs == UserDefFunction) {
        cRL = HALF_F*(csIsoArr(k,j,i)+csIsoArr(k-koffset,j-joffset,i-ioffset));
      } else {
        cRL = csIso;
      }
#endif
      cmax = FMAX(FABS(vRL[Xn]+cRL),FABS(vRL[Xn]-cRL));


      // 3-- Compute the conservative variables
      K_PrimToCons(uL, vL, gamma_m1);
      K_PrimToCons(uR, vR, gamma_m1);

      // 4-- Compute the left and right fluxes
      K_Flux(fluxL, vL, uL, cRL*cRL, Xn);
      K_Flux(fluxR, vR, uR, cRL*cRL, Xn);

      // 5-- Compute the flux from the left and right states
#pragma unroll
      for(int nv = 0 ; nv < NVAR; nv++) {
        Flux(nv,k,j,i) = HALF_F*(fluxL[nv]+fluxR[nv] - cmax*(uR[nv]-uL[nv]));
      }

      //6-- Compute maximum wave speed for this sweep
      cMax(k,j,i) = cmax;
    }
  );

  idfx::popRegion();
}

#endif // HYDRO_HDSOLVERS_TVDLFHD_HPP_
