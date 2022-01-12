// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_FLUXHD_HPP_
#define HYDRO_FLUXHD_HPP_
#include "idefix.hpp"
#include "hydro.hpp"


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

#endif //HYDRO_FLUXHD_HPP_
