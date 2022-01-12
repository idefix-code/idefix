// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_FLUXMHD_HPP_
#define HYDRO_FLUXMHD_HPP_

#include "idefix.hpp"
#include "hydro.hpp"

/********************************************************************************************
 * @fn void K_Flux(real F[], real V[], real U[], real Cs2Iso,
 *                                  const int Xn, const int Xt, const int Xb,
 *                                  const int BXn, const int BXt, const int BXb)
 * @param F[]   Array of flux variables (output)
 * @param V[]   Array of primitive variabless (input)
 * @param U[]   Array of conservative variables (input)
 * @param cs2Iso Isothermal sound speed (only used when ISOTHERMAL is defined)
 * @param Xn    Index of the normal velocity component
 * @param Xt    Index of the tangential velocity component
 * @param Xb    Indefx of the second tangential velocity component
 *
 *  This routine computes the MHD out of V and U variables and stores it in F
 ********************************************************************************************/
KOKKOS_INLINE_FUNCTION void K_Flux(real F[], real V[], real U[], real Cs2Iso,
                                   ARG_EXPAND(const int Xn, const int Xt, const int Xb),
                                   ARG_EXPAND(const int BXn, const int BXt, const int BXb)) {
  F[RHO] = U[Xn];
  EXPAND( F[MX1] = U[MX1]*V[Xn] - V[BXn]*V[BX1]; ,
          F[MX2] = U[MX2]*V[Xn] - V[BXn]*V[BX2]; ,
          F[MX3] = U[MX3]*V[Xn] - V[BXn]*V[BX3];)

  EXPAND(F[BXn] = ZERO_F;                             ,
          F[BXt] = V[Xn]*V[BXt] - V[BXn]*V[Xt];   ,
          F[BXb] = V[Xn]*V[BXb] - V[BXn]*V[Xb]; )

  real Bmag2 = EXPAND(V[BX1]*V[BX1] , + V[BX2]*V[BX2], + V[BX3]*V[BX3]);

#if HAVE_ENERGY
  real ptot  = V[PRS] + HALF_F*Bmag2;

#elif defined(ISOTHERMAL)
  real ptot  = Cs2Iso * V[RHO] + HALF_F*Bmag2;

#else
  #error "K_Flux not defined for this EOS!"
#endif

#if HAVE_ENERGY
  F[ENG]   = (U[ENG] + ptot)*V[Xn] - V[BXn] * (EXPAND( V[VX1]*V[BX1]   ,
                                                      + V[VX2]*V[BX2]  ,
                                                      + V[VX3]*V[BX3]  ));
#endif

  // Add back pressure in the flux (not included in original PLUTO implementation)
  F[Xn]   += ptot;
}

#endif // HYDRO_FLUXMHD_HPP_
