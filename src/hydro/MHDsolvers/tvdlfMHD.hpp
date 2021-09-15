// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_MHDSOLVERS_TVDLFMHD_HPP_
#define HYDRO_MHDSOLVERS_TVDLFMHD_HPP_

#include "../idefix.hpp"
#include "extrapolatePrimVar.hpp"
#include "fluxMHD.hpp"
#include "convertConsToPrimMHD.hpp"
#include "storeFlux.hpp"

// Compute Riemann fluxes from states using TVDLF solver
template<const int DIR>
void Hydro::TvdlfMHD() {
  idfx::pushRegion("Hydro::TVDLF_MHD");

  int ioffset,joffset,koffset;
  int iextend, jextend,kextend;

  ioffset=joffset=koffset=0;
  // extension in perp to the direction of integration, as required by CT.
  iextend=jextend=kextend=0;

  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Vs = this->Vs;
  IdefixArray4D<real> Flux = this->FluxRiemann;
  IdefixArray3D<real> cMax = this->cMax;
  IdefixArray3D<real> csIsoArr = this->isoSoundSpeedArray;

  // Required for high order interpolations
  IdefixArray1D<real> dx = this->data->dx[DIR];

  // References to required emf components
  IdefixArray3D<real> Eb;
  IdefixArray3D<real> Et;

#if EMF_AVERAGE == UCT_CONTACT
  IdefixArray3D<int> SV;
#elif EMF_AVERAGE == UCT_HLLD || EMF_AVERAGE == UCT_HLL2
  IdefixArray3D<real> aL;
  IdefixArray3D<real> aR;
  IdefixArray3D<real> dL;
  IdefixArray3D<real> dR;
#endif

  real gamma = this->gamma;
  real gamma_m1=gamma-ONE_F;
  real csIso = this->isoSoundSpeed;
  HydroModuleStatus haveIsoCs = this->haveIsoSoundSpeed;

  // Define normal, tangent and bi-tanget indices
  // st and sb will be useful only when Hall is included
  real st,sb;

  switch(DIR) {
    case(IDIR):
      ioffset = 1;
      D_EXPAND(
                st = -ONE_F;  ,
                jextend = 1;  ,
                kextend = 1;
                sb = +ONE_F;  )

      Et = this->emf.ezi;
      Eb = this->emf.eyi;
#if EMF_AVERAGE == UCT_CONTACT
      SV = this->emf.svx;
#elif EMF_AVERAGE == UCT_HLLD || EMF_AVERAGE == UCT_HLL
      aL = this->emf.axL;
      aR = this->emf.axR;

      dL = this->emf.dxL;
      dR = this->emf.dxR;
#endif
      break;
#if DIMENSIONS >= 2
    case(JDIR):
      joffset=1;
      D_EXPAND(
                iextend = 1;
                st = +ONE_F;  ,
                              ,
                kextend = 1;
                sb = -ONE_F;  )

      Et = this->emf.ezj;
      Eb = this->emf.exj;
  #if EMF_AVERAGE == UCT_CONTACT
      SV = this->emf.svy;
  #elif EMF_AVERAGE == UCT_HLLD || EMF_AVERAGE == UCT_HLL
      aL = this->emf.ayL;
      aR = this->emf.ayR;

      dL = this->emf.dyL;
      dR = this->emf.dyR;
  #endif
      break;
#endif
#if DIMENSIONS == 3
    case(KDIR):
      koffset=1;
      D_EXPAND(
                iextend = 1;
                st = -ONE_F;  ,
                jextend = 1;  ,
                sb = +ONE_F;  )

      Et = this->emf.eyk;
      Eb = this->emf.exk;
  #if EMF_AVERAGE == UCT_CONTACT
      SV = this->emf.svz;
  #elif EMF_AVERAGE == UCT_HLLD || EMF_AVERAGE == UCT_HLL
      aL = this->emf.azL;
      aR = this->emf.azR;

      dL = this->emf.dzL;
      dR = this->emf.dzR;
  #endif
      break;
#endif
    default:
      IDEFIX_ERROR("Wrong direction");
  }

  idefix_for("CalcRiemannFlux",
             data->beg[KDIR]-kextend,data->end[KDIR]+koffset+kextend,
             data->beg[JDIR]-jextend,data->end[JDIR]+joffset+jextend,
             data->beg[IDIR]-iextend,data->end[IDIR]+ioffset+iextend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      // Init the directions (should be in the kernel for proper optimisation by the compilers)
      EXPAND( const int Xn = DIR+MX1;                    ,
              const int Xt = (DIR == IDIR ? MX2 : MX1);  ,
              const int Xb = (DIR == KDIR ? MX2 : MX3);  )

      EXPAND( const int BXn = DIR+BX1;                    ,
              const int BXt = (DIR == IDIR ? BX2 : BX1);  ,
              const int BXb = (DIR == KDIR ? BX2 : BX3);   )
      // Primitive variables
      real vL[NVAR];
      real vR[NVAR];
      real v[NVAR];

      real uL[NVAR];
      real uR[NVAR];

      real fluxL[NVAR];
      real fluxR[NVAR];

      // Load primitive variables
      K_ExtrapolatePrimVar<DIR>(i, j, k, Vc, Vs, dx, vL, vR);
#pragma unroll
      for(int nv = 0 ; nv < NVAR; nv++) {
        v[nv] = HALF_F*(vL[nv] + vR[nv]);
      }

      // Get the wave speed
      // Signal speeds
      real cRL, cmax, c2Iso;
      real gpr, Bt2, B2;

      // Init c2Isothermal (used only when isothermal approx is set)
      c2Iso = ZERO_F;

#if HAVE_ENERGY
      gpr=gamma*v[PRS];
#else
      if(haveIsoCs == UserDefFunction) {
        c2Iso = HALF_F*(csIsoArr(k,j,i)+csIsoArr(k-koffset,j-joffset,i-ioffset));
        c2Iso = c2Iso*c2Iso;
      } else {
        c2Iso = csIso*csIso;
      }

      gpr = c2Iso*v[RHO];
#endif
      Bt2=EXPAND( ZERO_F           ,
                  + v[BXt]*v[BXt]  ,
                  + v[BXb]*v[BXb]  );

      B2=Bt2 + v[BXn]*v[BXn];

      cRL = gpr - B2;
      cRL = cRL + B2 + std::sqrt(cRL*cRL + FOUR_F*gpr*Bt2);
      cRL = std::sqrt(HALF_F * cRL/v[RHO]);

      cmax = std::fmax(std::fabs(v[Xn]+cRL),FABS(v[Xn]-cRL));

      real sl, sr;
      sl = -cmax;
      sr = cmax;


      // 2-- Compute the conservative variables
      K_PrimToCons(uL, vL, gamma_m1);
      K_PrimToCons(uR, vR, gamma_m1);

      // 3-- Compute the left and right fluxes
      K_Flux(fluxL, vL, uL, c2Iso, ARG_EXPAND(Xn, Xt, Xb), ARG_EXPAND(BXn, BXt, BXb));
      K_Flux(fluxR, vR, uR, c2Iso, ARG_EXPAND(Xn, Xt, Xb), ARG_EXPAND(BXn, BXt, BXb));


      // 5-- Compute the flux from the left and right states
#pragma unroll
      for(int nv = 0 ; nv < NVAR; nv++) {
        Flux(nv,k,j,i) = HALF_F*(fluxL[nv] + fluxR[nv] + cmax*(uL[nv] - uR[nv]));
      }

      //6-- Compute maximum wave speed for this sweep
      cMax(k,j,i) = cmax;

      // 7-- Store the flux in the emf components
      #if EMF_AVERAGE == ARITHMETIC || EMF_AVERAGE == UCT0
        K_StoreEMF<DIR>(i,j,k,st,sb, Flux, Et, Eb);
      #elif EMF_AVERAGE == UCT_CONTACT
        K_StoreContact<DIR>(i,j,k,st,sb,Flux,Et,Eb,SV);
      #elif EMF_AVERAGE == UCT_HLL2
        K_StoreHLL<DIR>(i,j,k,st,sb,sl,sr,vL,vR,Et,Eb,aL,aR,dL,dR);
      #elif EMF_AVERAGE == UCT_HLLD
        // Compute magnetic pressure
        real ptR, ptL;

        #if HAVE_ENERGY
          ptL  = vL[PRS] + HALF_F* ( EXPAND(vL[BX1]*vL[BX1]     ,
                                            + vL[BX2]*vL[BX2]   ,
                                            + vL[BX3]*vL[BX3])  );
          ptR  = vR[PRS] + HALF_F* ( EXPAND(vR[BX1]*vR[BX1]     ,
                                            + vR[BX2]*vR[BX2]   ,
                                            + vR[BX3]*vR[BX3])  );
        #else
          ptL  = c2Iso*vL[RHO] + HALF_F* (EXPAND(vL[BX1]*vL[BX1]     ,
                                                + vL[BX2]*vL[BX2]   ,
                                                + vL[BX3]*vL[BX3])  );
          ptR  = c2Iso*vR[RHO] + HALF_F* (EXPAND(vR[BX1]*vR[BX1]     ,
                                                + vR[BX2]*vR[BX2]   ,
                                                + vR[BX3]*vR[BX3])  );
        #endif
        K_StoreHLLD<DIR>(i,j,k,st,sb,sl,sr,ptL,ptR,vL,vR,uL,uR,Et,Eb,aL,aR,dL,dR);
      #else
        #error "Unknown EMF_AVERAGE scheme"
      #endif
  });

  idfx::popRegion();
}

#endif // HYDRO_MHDSOLVERS_TVDLFMHD_HPP_
