// ********************************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ********************************************************************************************************

#pragma once
#include "../idefix.hpp"

// Local Kokkos Inline functions

KOKKOS_INLINE_FUNCTION void K_Flux(real F[], real V[], real U[], real C2Iso, 
                                    const int Xn, const int Xt, const int Xb,
                                    const int BXn, const int BXt, const int BXb) {

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
    real ptot  = C2Iso * V[RHO] + HALF_F*Bmag2;
    #else
    #error "K_Flux not defined for this EOS!"
    #endif

    #if HAVE_ENERGY
    F[ENG]   = (U[ENG] + ptot)*V[Xn] - V[BXn] * (EXPAND(V[VX1]*V[BX1] , + V[VX2]*V[BX2], + V[VX3]*V[BX3]));
    #endif

    // Add back pressure in the flux (not included in original PLUTO implementation)
    F[Xn]   += ptot;
} 

KOKKOS_INLINE_FUNCTION void K_ConsToPrim(real Vc[], real Uc[], real gamma_m1) {
    

    Vc[RHO] = Uc[RHO];

    EXPAND( Vc[VX1] = Uc[MX1]/Uc[RHO]; ,
            Vc[VX2] = Uc[MX2]/Uc[RHO]; ,
            Vc[VX3] = Uc[MX3]/Uc[RHO];)

    EXPAND( Vc[BX1] = Uc[BX1];  ,
            Vc[BX2] = Uc[BX2];  ,
            Vc[BX3] = Uc[BX3];  )
         

#if HAVE_ENERGY
    real kin,mag;
    kin = HALF_F / Uc[RHO] * (EXPAND(    Uc[MX1]*Uc[MX1] , 
                                    + Uc[MX2]*Uc[MX2] ,
                                    + Uc[MX3]*Uc[MX3] ));

    mag = HALF_F * (EXPAND(    Uc[BX1]*Uc[BX1] , 
                                    + Uc[BX2]*Uc[BX2] ,
                                    + Uc[BX3]*Uc[BX3]));


    Vc[PRS] = gamma_m1 * (Uc[ENG] - kin - mag);
    
    // Check pressure positivity
    if(Vc[PRS]<= ZERO_F) {
      #ifdef SMALL_PRESSURE_TEMPERATURE
        Vc[PRS] = SMALL_PRESSURE_TEMPERATURE*Vc_RHO;
        #else
        Vc[PRS] = SMALL_PRESSURE_FIX;
        #endif
        
        Uc[ENG] = Vc[PRS]/gamma_m1+kin+mag;
    }
#endif  // Have_energy
}

KOKKOS_INLINE_FUNCTION void K_PrimToCons(real Uc[], real Vc[], real gamma_m1) {

    Uc[RHO] = Vc[RHO];

    EXPAND( Uc[MX1] = Vc[VX1]*Vc[RHO]; ,
            Uc[MX2] = Vc[VX2]*Vc[RHO]; ,
            Uc[MX3] = Vc[VX3]*Vc[RHO];)


    EXPAND( Uc[BX1] = Vc[BX1];  ,
            Uc[BX2] = Vc[BX2];  ,
            Uc[BX3] = Vc[BX3];  )

#if HAVE_ENERGY

    Uc[ENG] = Vc[PRS] / gamma_m1 
                + HALF_F * Vc[RHO] * (EXPAND(  Vc[VX1]*Vc[VX1] , 
                                         + Vc[VX2]*Vc[VX2] ,
                                         + Vc[VX3]*Vc[VX3] ))
                + HALF_F * (EXPAND(    Uc[BX1]*Uc[BX1] , 
                                     + Uc[BX2]*Uc[BX2] ,
                                     + Uc[BX3]*Uc[BX3])); 
#endif  // Have_energy

}

#include "MHDsolvers/tvdlfMHD.hpp"
#include "MHDsolvers/hllMHD.hpp"
#include "MHDsolvers/hlldMHD.hpp"
#include "MHDsolvers/roeMHD.hpp"
