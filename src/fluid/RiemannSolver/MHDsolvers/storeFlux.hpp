// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_RIEMANNSOLVER_MHDSOLVERS_STOREFLUX_HPP_
#define FLUID_RIEMANNSOLVER_MHDSOLVERS_STOREFLUX_HPP_

#include "idefix.hpp"
#include "fluid.hpp"

template <const int DIR>
KOKKOS_FORCEINLINE_FUNCTION void K_StoreEMF( const int i, const int j, const int k,
                                        const real st, const real sb,
                                        const IdefixArray4D<real> &Flux,
                                        const IdefixArray3D<real> &Et,
                                        const IdefixArray3D<real> &Eb ) {
  EXPAND(                                           ,
                         constexpr int BXt = (DIR == IDIR ? BX2 : BX1);  ,
        [[maybe_unused]] constexpr int BXb = (DIR == KDIR ? BX2 : BX3);   )

  #if COMPONENTS > 1
    D_EXPAND( Et(k,j,i) = st*Flux(BXt,k,j,i);  ,
                                                  ,
              Eb(k,j,i) = sb*Flux(BXb,k,j,i);  )
  #endif
}

template <const int DIR>
KOKKOS_FORCEINLINE_FUNCTION void K_StoreContact( const int i, const int j, const int k,
                                        const real st, const real sb,
                                        const IdefixArray4D<real> &Flux,
                                        const IdefixArray3D<real> &Et,
                                        const IdefixArray3D<real> &Eb,
                                        const IdefixArray3D<real> &SV) {
  K_StoreEMF<DIR>(i,j,k,st,sb,Flux,Et,Eb);
  real s = HALF_F;
  if (Flux(RHO,k,j,i) >  eps_UCT_CONTACT) s =  ONE_F;
  if (Flux(RHO,k,j,i) < -eps_UCT_CONTACT) s = ZERO_F;

  SV(k,j,i) = s;
}

template <const int DIR>
KOKKOS_FORCEINLINE_FUNCTION void K_StoreHLL( const int i, const int j, const int k,
                                        const real st, const real sb,
                                        const real sl, const real sr,
                                        real vL[], real vR[],
                                        const IdefixArray3D<real> &Et,
                                        const IdefixArray3D<real> &Eb,
                                        const IdefixArray3D<real> &aL,
                                        const IdefixArray3D<real> &aR,
                                        const IdefixArray3D<real> &dL,
                                        const IdefixArray3D<real> &dR) {
  EXPAND(                                          ,
        constexpr int Xt = (DIR == IDIR ? MX2 : MX1);  ,
        constexpr int Xb = (DIR == KDIR ? MX2 : MX3);  )

  real ar = std::fmax(ZERO_F, sr);
  real al = std::fmin(ZERO_F, sl);
  real scrh = ONE_F/(ar - al);

  #if COMPONENTS > 1
  EXPAND( Et(k,j,i) = -st*(ar*vL[Xt] - al*vR[Xt])*scrh;  ,
                                                        ,
          Eb(k,j,i) = -sb*(ar*vL[Xb] - al*vR[Xb])*scrh;  );
  #endif

  aL(k,j,i) =  ar*scrh;
  aR(k,j,i) = -al*scrh;
  dR(k,j,i) = -al*ar*scrh;
  dL(k,j,i) =  dR(k,j,i);
}

template <const int DIR>
KOKKOS_FORCEINLINE_FUNCTION void K_StoreHLLD( const int i, const int j, const int k,
                                        const real st, const real sb,
                                        const real c2Iso,
                                        const real sl, const real sr,
                                        real vL[], real vR[],
                                        real uL[], real uR[],
                                        const IdefixArray3D<real> &Et,
                                        const IdefixArray3D<real> &Eb,
                                        const IdefixArray3D<real> &aL,
                                        const IdefixArray3D<real> &aR,
                                        const IdefixArray3D<real> &dL,
                                        const IdefixArray3D<real> &dR) {
  EXPAND( const int Xn = DIR+MX1;                    ,
        const int Xt = (DIR == IDIR ? MX2 : MX1);  ,
        const int Xb = (DIR == KDIR ? MX2 : MX3);  )
  // Compute magnetic pressure
  [[maybe_unused]] real ptR, ptL;

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

  int revert_to_hll = 0;

  const int BXn = DIR+BX1;

  real Bn = (sr*vR[BXn] - sl*vL[BXn])/(sr - sl);

  real chiL, chiR, nuLR, nuL, nuR, SaL, SaR;
  real Sc;
  real eps = 1.e-12*(fabs(sl) + fabs(sr));
  real duL  = sl - vL[Xn];
  real duR  = sr - vR[Xn];

#if HAVE_ENERGY
  // Recompute speeds
  real sqrL, sqrR, usLRHO, usRRHO;

  real scrh  = ONE_F/(duR*uR[RHO] - duL*uL[RHO]);
  Sc = (duR*uR[Xn] - duL*uL[Xn] - ptR + ptL)*scrh;

  usLRHO = uL[RHO]*duL/(sl - Sc);
  usRRHO = uR[RHO]*duR/(sr - Sc);

  sqrL = sqrt(usLRHO);
  sqrR = sqrt(usRRHO);

  SaL = Sc - fabs(Bn)/sqrL;
  SaR = Sc + fabs(Bn)/sqrR;

  if((SaL - sl) < 1e-4*(Sc-sl)) revert_to_hll = 1;
  if((SaR - sr) > -1e-4*(sr-Sc)) revert_to_hll = 1;

  chiL  = (vL[Xn] - Sc)*(sl - Sc)/(SaL + sl - TWO_F*Sc);
  chiR  = (vR[Xn] - Sc)*(sr - Sc)/(SaR + sr - TWO_F*Sc);
#else
  real scrh    = ONE_F/(sr - sl);
  real rho_h   = (uR[RHO]*duR - uL[RHO]*duL)*scrh;
  Sc = (sl*uR[RHO]*duR - sr*uL[RHO]*duL)*scrh/rho_h;
  // Recompute speeds
  real sqrho_h = sqrt(rho_h);
  SaL = Sc - fabs(Bn)/sqrho_h;
  SaR = Sc + fabs(Bn)/sqrho_h;

  if((SaL - sl) < 1e-4*(Sc-sl)) revert_to_hll = 1;
  if((SaR - sr) > -1e-4*(sr-Sc)) revert_to_hll = 1;

  chiL  = (vL[Xn] - Sc)*(sl - Sc)/(SaL + sl - TWO_F*Sc);
  chiR  = (vR[Xn] - Sc)*(sr - Sc)/(SaR + sr - TWO_F*Sc);
#endif

  nuL  = (SaL + sl)/(fabs(SaL) + fabs(sl));
  nuR  = (SaR + sr)/(fabs(SaR) + fabs(sr));
  nuLR = (SaL + SaR)/(fabs(SaL) + fabs(SaR));

  if (fabs(SaR - SaL) > 1.e-9*fabs(sr-sl)) {
    dL(k,j,i) =   HALF_F*(chiL*nuL - chiL*nuLR)
                    + HALF_F*(fabs(SaL) - nuLR*SaL);

    dR(k,j,i) =   HALF_F*(chiR*nuR - chiR*nuLR)
                    + HALF_F*(fabs(SaR) - nuLR*SaR);
    aL(k,j,i) = HALF_F*(ONE_F + nuLR);
    aR(k,j,i) = HALF_F*(ONE_F - nuLR);

  } else {   // HLLC, degenerate limit Bx -> 0
    dL(k,j,i) = HALF_F*chiL*nuL + HALF_F*fabs(SaL);
    dR(k,j,i) = HALF_F*chiR*nuR + HALF_F*fabs(SaR);

    aL(k,j,i) = HALF_F;
    aR(k,j,i) = HALF_F;
  }

  real ar = std::fmax(ZERO_F, sr);
  real al = std::fmin(ZERO_F, sl);
  scrh = ONE_F/(ar - al);

  // HLL diffusion coefficients
  if (revert_to_hll) {
    aL(k,j,i) =  ar*scrh;
    aR(k,j,i) = -al*scrh;
    dR(k,j,i) = -al*ar*scrh;
    dL(k,j,i) =  dR(k,j,i);
  }

  // Lax-Friedrichs diffusion coefficients
  if(0) {
    real lambda = std::fmax(sr,sl);
    aL(k,j,i) = HALF_F;
    aR(k,j,i) = HALF_F;
    dR(k,j,i) = HALF_F*lambda;
    dL(k,j,i) = HALF_F*lambda;
  }

  #if COMPONENTS > 1
  EXPAND( Et(k,j,i) = -st*(ar*vL[Xt] - al*vR[Xt])*scrh;  ,
                                                          ,
          Eb(k,j,i) = -sb*(ar*vL[Xb] - al*vR[Xb])*scrh;  );
  #endif
}

#endif //FLUID_RIEMANNSOLVER_MHDSOLVERS_STOREFLUX_HPP_
