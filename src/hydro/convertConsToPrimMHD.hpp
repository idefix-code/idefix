// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_CONVERTCONSTOPRIMMHD_HPP_
#define HYDRO_CONVERTCONSTOPRIMMHD_HPP_

#include "idefix.hpp"
#include "hydro.hpp"


KOKKOS_INLINE_FUNCTION void K_ConsToPrim(real Vc[], real Uc[], real gamma_m1) {
  Vc[RHO] = Uc[RHO];

  EXPAND( Vc[VX1] = Uc[MX1]/Uc[RHO];  ,
          Vc[VX2] = Uc[MX2]/Uc[RHO];  ,
          Vc[VX3] = Uc[MX3]/Uc[RHO];  )

  EXPAND( Vc[BX1] = Uc[BX1];  ,
          Vc[BX2] = Uc[BX2];  ,
          Vc[BX3] = Uc[BX3];  )


#if HAVE_ENERGY
  real kin,mag;
  kin = HALF_F / Uc[RHO] * (EXPAND( Uc[MX1]*Uc[MX1]   ,
                                   + Uc[MX2]*Uc[MX2]  ,
                                   + Uc[MX3]*Uc[MX3]  ));

  mag = HALF_F * (EXPAND( Uc[BX1]*Uc[BX1]   ,
                         + Uc[BX2]*Uc[BX2]  ,
                         + Uc[BX3]*Uc[BX3]  ));


  Vc[PRS] = gamma_m1 * (Uc[ENG] - kin - mag);

  // Check pressure positivity
  if(Vc[PRS]<= ZERO_F) {
  #ifdef SMALL_PRESSURE_TEMPERATURE
    Vc[PRS] = SMALL_PRESSURE_TEMPERATURE*Vc[RHO];
  #else
    Vc[PRS] = SMALL_PRESSURE_FIX;
  #endif

    Uc[ENG] = Vc[PRS]/gamma_m1+kin+mag;
  }
#endif  // Have_energy
}

KOKKOS_INLINE_FUNCTION void K_PrimToCons(real Uc[], real Vc[], real gamma_m1) {
  Uc[RHO] = Vc[RHO];

  EXPAND( Uc[MX1] = Vc[VX1]*Vc[RHO];  ,
          Uc[MX2] = Vc[VX2]*Vc[RHO];  ,
          Uc[MX3] = Vc[VX3]*Vc[RHO];  )


  EXPAND( Uc[BX1] = Vc[BX1];  ,
          Uc[BX2] = Vc[BX2];  ,
          Uc[BX3] = Vc[BX3];  )

#if HAVE_ENERGY

  Uc[ENG] = Vc[PRS] / gamma_m1
              + HALF_F * Vc[RHO] * (EXPAND( Vc[VX1]*Vc[VX1]  ,
                                          + Vc[VX2]*Vc[VX2]  ,
                                          + Vc[VX3]*Vc[VX3]  ))
              + HALF_F * (EXPAND( Uc[BX1]*Uc[BX1]  ,
                                + Uc[BX2]*Uc[BX2]  ,
                                + Uc[BX3]*Uc[BX3]  ));
#endif  // Have_energy
}




#endif // HYDRO_CONVERTCONSTOPRIMMHD_HPP_
