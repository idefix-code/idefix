// ***********************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************************

#ifndef HYDRO_MHDSOLVERS_TVDLFMHD_HPP_
#define HYDRO_MHDSOLVERS_TVDLFMHD_HPP_

#include "../idefix.hpp"
#include "solversMHD.hpp"

// Compute Riemann fluxes from states using TVDLF solver
template<const int DIR, ARG_EXPAND(const int Xn, const int Xt, const int Xb),
                        ARG_EXPAND(const int BXn, const int BXt, const int BXb)>
void Hydro::TvdlfMHD() {
  idfx::pushRegion("Hydro::TVDLF_MHD");

  int ioffset,joffset,koffset;
  int iextend, jextend,kextend;

  ioffset=joffset=koffset=0;
  // extension in perp to the direction of integration, as required by CT.
  iextend=jextend=kextend=0;

  IdefixArray4D<real> PrimL = this->PrimL;
  IdefixArray4D<real> PrimR = this->PrimR;
  IdefixArray4D<real> Flux = this->FluxRiemann;
  IdefixArray3D<real> cMax = this->cMax;

  // References to required emf components
  IdefixArray3D<real> Eb;
  IdefixArray3D<real> Et;

  IdefixArray3D<int> SV;

  real gamma = this->gamma;
  real gamma_m1=gamma-ONE_F;
  real C2Iso = this->C2Iso;

  // Define normal, tangent and bi-tanget indices
  // st and sb will be useful only when Hall is included
  D_EXPAND( real st;  ,
                      ,
            real sb;  )

  switch(DIR) {
    case(IDIR):
      ioffset = 1;
      D_EXPAND(               ,
                jextend = 1;  ,
                kextend = 1;  )

      Et = this->emf.ezi;
      Eb = this->emf.eyi;
      SV = this->emf.svx;

      D_EXPAND( st = -1.0;  ,
                            ,
                sb = +1.0;  )
      break;
    case(JDIR):
      joffset=1;
      D_EXPAND( iextend = 1;  ,
                              ,
                kextend = 1;  )

      Et = this->emf.ezj;
      Eb = this->emf.exj;
      SV = this->emf.svy;

      D_EXPAND( st = +1.0;  ,
                            ,
                sb = -1.0;  )
      break;
    case(KDIR):
      koffset=1;
      D_EXPAND( iextend = 1;  ,
                jextend = 1;  ,
                              )

      Et = this->emf.eyk;
      Eb = this->emf.exk;
      SV = this->emf.svz;

      D_EXPAND( st = -1.0;  ,
                            ,
                sb = +1.0;  )
      break;
    default:
      IDEFIX_ERROR("Wrong direction");
  }

  idefix_for("CalcRiemannFlux",
             data->beg[KDIR]-kextend,data->end[KDIR]+koffset+kextend,
             data->beg[JDIR]-jextend,data->end[JDIR]+joffset+jextend,
             data->beg[IDIR]-iextend,data->end[IDIR]+ioffset+iextend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      // Primitive variables
      real vL[NVAR];
      real vR[NVAR];
      real v[NVAR];

      real uL[NVAR];
      real uR[NVAR];

      real fluxL[NVAR];
      real fluxR[NVAR];

      // Load primitive variables
#pragma unroll
      for(int nv = 0 ; nv < NVAR; nv++) {
        vL[nv] = PrimL(nv,k,j,i);
        vR[nv] = PrimR(nv,k,j,i);
        v[nv] = HALF_F*(vL[nv] + vR[nv]);
      }

      // Get the wave speed
      // Signal speeds
      real cRL, cmax;
      real gpr, Bt2, B2;

#if HAVE_ENERGY
      gpr=gamma*v[PRS];
#else
      gpr=C2Iso*v[RHO];
#endif
      Bt2=EXPAND( ZERO_F           ,
                  + v[BXt]*v[BXt]  ,
                  + v[BXb]*v[BXb]  );

      B2=Bt2 + v[BXn]*v[BXn];

      cRL = gpr - B2;
      cRL = cRL + B2 + std::sqrt(cRL*cRL + FOUR_F*gpr*Bt2);
      cRL = std::sqrt(HALF_F * cRL/v[RHO]);

      cmax = FMAX(FABS(v[Xn]+cRL),FABS(v[Xn]-cRL));


      // 2-- Compute the conservative variables
      K_PrimToCons(uL, vL, gamma_m1);
      K_PrimToCons(uR, vR, gamma_m1);

      // 3-- Compute the left and right fluxes
      K_Flux(fluxL, vL, uL, C2Iso, ARG_EXPAND(Xn, Xt, Xb), ARG_EXPAND(BXn, BXt, BXb));
      K_Flux(fluxR, vR, uR, C2Iso, ARG_EXPAND(Xn, Xt, Xb), ARG_EXPAND(BXn, BXt, BXb));


      // 5-- Compute the flux from the left and right states
#pragma unroll
      for(int nv = 0 ; nv < NVAR; nv++) {
        Flux(nv,k,j,i) = HALF_F*(fluxL[nv] + fluxR[nv] + cmax*(uL[nv] - uR[nv]));
      }

      //6-- Compute maximum wave speed for this sweep
      cMax(k,j,i) = cmax;

      // 7-- Store the flux in the emf components
      D_EXPAND( Et(k,j,i) = st*Flux(BXt,k,j,i);  ,
                                                 ,
                Eb(k,j,i) = sb*Flux(BXb,k,j,i);  )

#if EMF_AVERAGE == UCT_CONTACT
      int s = 0;
      if (Flux(RHO,k,j,i) >  eps_UCT_CONTACT) s =  1;
      if (Flux(RHO,k,j,i) < -eps_UCT_CONTACT) s = -1;

      SV(k,j,i) = s;
#endif
  });

  idfx::popRegion();
}

#endif // HYDRO_MHDSOLVERS_TVDLFMHD_HPP_
