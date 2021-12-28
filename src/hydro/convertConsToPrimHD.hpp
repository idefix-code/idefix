// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_CONVERTCONSTOPRIMHD_HPP_
#define HYDRO_CONVERTCONSTOPRIMHD_HPP_

#include "idefix.hpp"
#include "hydro.hpp"

KOKKOS_INLINE_FUNCTION void K_ConsToPrim(
                                  real Vc[],
                                  real Uc[],
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

  // Check pressure positivity
  if(Vc[PRS]<= ZERO_F) {
  #ifdef SMALL_PRESSURE_TEMPERATURE
    Vc[PRS] = SMALL_PRESSURE_TEMPERATURE*Vc[RHO];
  #else
    Vc[PRS] = SMALL_PRESSURE_FIX;
  #endif
    Uc[ENG] = Vc[PRS]/gamma_m1+kin;
  }
#endif  // Have_energy
}

KOKKOS_INLINE_FUNCTION void K_PrimToCons(real *KOKKOS_RESTRICT Uc,
                                    const real *KOKKOS_RESTRICT Vc,
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

#endif // HYDRO_CONVERTCONSTOPRIMHD_HPP_
