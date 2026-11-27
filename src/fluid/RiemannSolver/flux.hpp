// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/fluid/RiemannSolver/flux.hpp
//
// Last modified : 10/2023
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2023)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_RIEMANNSOLVER_FLUX_HPP_
#define FLUID_RIEMANNSOLVER_FLUX_HPP_
#include "idefix.hpp"
#include "fluid.hpp"


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
template<typename Phys, int DIR>
KOKKOS_INLINE_FUNCTION void K_Flux(real *KOKKOS_RESTRICT F, const real *KOKKOS_RESTRICT V,
                                   const real *KOKKOS_RESTRICT U, real Cs2Iso) {
  constexpr int Xn = DIR+MX1;
  [[maybe_unused]] constexpr int BXn = DIR+BX1;

  // Mass flux (common to all physics)
  F[RHO] = U[Xn];

  // Momentum flux
  if constexpr(Phys::mhd) {
    EXPAND( F[MX1] = U[MX1]*V[Xn] - V[BXn]*V[BX1]; ,
            F[MX2] = U[MX2]*V[Xn] - V[BXn]*V[BX2]; ,
            F[MX3] = U[MX3]*V[Xn] - V[BXn]*V[BX3];)
  } else {
    EXPAND ( F[MX1] = U[MX1]*V[Xn];  ,
             F[MX2] = U[MX2]*V[Xn];  ,
             F[MX3] = U[MX3]*V[Xn];  )
  }

  // Magnetic flux
  if constexpr(Phys::mhd) {
    EXPAND(                                                 ,
            constexpr int Xt = (DIR == IDIR ? MX2 : MX1);  ,
            constexpr int Xb = (DIR == KDIR ? MX2 : MX3);  )

    EXPAND(                                                 ,
            constexpr int BXt = (DIR == IDIR ? BX2 : BX1);  ,
            constexpr int BXb = (DIR == KDIR ? BX2 : BX3);   )

    EXPAND(F[BXn] = ZERO_F;                             ,
            F[BXt] = V[Xn]*V[BXt] - V[BXn]*V[Xt];   ,
            F[BXb] = V[Xn]*V[BXb] - V[BXn]*V[Xb]; )
  }

  if constexpr(Phys::pressure || Phys::isothermal) {
    // Pressure-related term
    real ptot;
    if constexpr(Phys::mhd) {
      ////////////////
      // MHD VERSION
      ///////////////
      real Bmag2 = EXPAND(V[BX1]*V[BX1] , + V[BX2]*V[BX2], + V[BX3]*V[BX3]);
      if constexpr(Phys::pressure) {
        ptot  = V[PRS] + HALF_F*Bmag2;
        // Energy flux
        F[ENG]   = (U[ENG] + ptot)*V[Xn] - V[BXn] * (EXPAND( V[VX1]*V[BX1]   ,
                                                      + V[VX2]*V[BX2]  ,
                                                      + V[VX3]*V[BX3]  ));

      } else if constexpr(Phys::isothermal) {
        ptot  = Cs2Iso * V[RHO] + HALF_F*Bmag2;
      }
    } else {
      ////////////////
      // Hydro VERSION
      ///////////////
      if constexpr(Phys::pressure) {
        ptot = V[PRS];
        F[ENG]  = (U[ENG] + ptot)*V[Xn];
      } else if constexpr(Phys::isothermal) {
        ptot  = Cs2Iso * V[RHO];
      }
    }
    // Add back pressure in the flux (not included in original PLUTO implementation)
    F[Xn]   += ptot;
  }
}


/********************************************************************************************
 * @fn void K_Flux(real F[], real V[], real U[], real reduced_c, real xi
 *                                  const int Xn, const int Xt, const int Xb)
 * @param F[]   Array of flux variables (output)
 * @param V[]   Array of primitive variabless (input)
 * @param U[]   Array of conservative variables (input)
 * @param reduced_c reduced speed of light for radiation
 * @param xi    Closure parameter for radiation pressure
 * @param Xn    Index of the normal velocity component
 *
 *  This routine computes the radiation out of V and U variables and stores it in F
 ********************************************************************************************/
template<typename Phys, int DIR>
KOKKOS_INLINE_FUNCTION void K_Flux(real *KOKKOS_RESTRICT F, const real *KOKKOS_RESTRICT V,
                                   const real *KOKKOS_RESTRICT U, real reduced_c, real xi) {
  constexpr int Xn = DIR+MX1;

  // Radiation energy flux
  F[ER] = U[Xn]*reduced_c;

  // Add momentum-like part of the radiation pressure tensor
  real Fnorm2 = EXPAND(V[FR1]*V[FR1] , + V[FR2]*V[FR2], + V[FR3]*V[FR3]);
  real inv_Fnorm2 = (Fnorm2 <= 1.e-100 ? 1.e-100 : ONE_F / Fnorm2);

  EXPAND ( F[FR1] = HALF_F*(3.*xi-1.)*V[ER]*V[FR1]*V[Xn]*inv_Fnorm2;  ,
           F[FR2] = HALF_F*(3.*xi-1.)*V[ER]*V[FR2]*V[Xn]*inv_Fnorm2;  ,
           F[FR3] = HALF_F*(3.*xi-1.)*V[ER]*V[FR3]*V[Xn]*inv_Fnorm2;  )

  // Add pressure-like part of the radiation pressure tensor
  F[Xn] += HALF_F*(1.-xi)*V[ER],

  EXPAND ( F[FR1] *= reduced_c; ,
           F[FR2] *= reduced_c; ,
           F[FR3] *= reduced_c; )
}

#endif //FLUID_RIEMANNSOLVER_FLUX_HPP_
