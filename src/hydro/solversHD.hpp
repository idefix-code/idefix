// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#pragma once
#include "idefix.hpp"


// Local Kokkos Inlined functions

/********************************************************************************************
 * @fn void K_Flux(real F[], real V[], real U[], real Cs2Iso,
 *                                  const int Xn, const int Xt, const int Xb,
 *                                  const int BXn, const int BXt, const int BXb)
 * @param F[]   Array of flux variables (output)
 * @param V[]   Array of primitive variabless (input)
 * @param U[]   Array of conservative variables (input)
 * @param cs2Iso Isothermal sound speed (only used when ISOTHERMAL is defined)
 * @param Xn    Index of the normal velocity component
 * 
 *  This routine computes the MHD out of V and U variables and stores it in F
 ********************************************************************************************/
KOKKOS_INLINE_FUNCTION void K_Flux(real *KOKKOS_RESTRICT F, const real *KOKKOS_RESTRICT V,
                                   const real *KOKKOS_RESTRICT U, real Cs2Iso, const int Xn) {
  F[RHO] = U[Xn];

  EXPAND( F[MX1] = U[MX1]*V[Xn];  ,
          F[MX2] = U[MX2]*V[Xn];  ,
          F[MX3] = U[MX3]*V[Xn];  )

#if HAVE_ENERGY
  F[ENG]  = (U[ENG] + V[PRS])*V[Xn];
  F[Xn]  += V[PRS];
#else
  // Add back pressure in the flux
  F[Xn]  += Cs2Iso * V[RHO];
#endif
}

KOKKOS_INLINE_FUNCTION void K_ConsToPrim(real *KOKKOS_RESTRICT Vc, const real *KOKKOS_RESTRICT Uc,
                                         real gamma_m1) {
  Vc[RHO] = Uc[RHO];

  EXPAND( Vc[VX1] = Uc[MX1]/Uc[RHO];  ,
          Vc[VX2] = Uc[MX2]/Uc[RHO];  ,
          Vc[VX3] = Uc[MX3]/Uc[RHO];  )

#if HAVE_ENERGY
  real kin;
  kin = HALF_F / Uc[RHO] * (EXPAND( Uc[MX1]*Uc[MX1]  ,
                                  + Uc[MX2]*Uc[MX2]  ,
                                  + Uc[MX3]*Uc[MX3]  ));

  Vc[PRS] = gamma_m1 * (Uc[ENG] - kin);
#endif  // Have_energy
}

KOKKOS_INLINE_FUNCTION void K_PrimToCons(real *KOKKOS_RESTRICT Uc, const real *KOKKOS_RESTRICT Vc,
                                         real gamma_m1) {
  Uc[RHO] = Vc[RHO];

  EXPAND( Uc[MX1] = Vc[VX1]*Vc[RHO];  ,
          Uc[MX2] = Vc[VX2]*Vc[RHO];  ,
          Uc[MX3] = Vc[VX3]*Vc[RHO];  )

#if HAVE_ENERGY

  Uc[ENG] = Vc[PRS] / gamma_m1
          + HALF_F * Vc[RHO] * (EXPAND(  Vc[VX1]*Vc[VX1]  ,
                                      + Vc[VX2]*Vc[VX2]   ,
                                      + Vc[VX3]*Vc[VX3]   ));
#endif  // Have_energy
}

#include "HDsolvers/tvdlfHD.hpp"
#include "HDsolvers/hllHD.hpp"
#include "HDsolvers/hllcHD.hpp"
#include "HDsolvers/roeHD.hpp"
