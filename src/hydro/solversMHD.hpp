// ********************************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ********************************************************************************************************

/********************************
 * Local Kokkos Inline function *
 * ******************************/
#pragma once
#include "../idefix.hpp"

void TvdlfMHD(DataBlock &, int, real, real);
void HllMHD(DataBlock &, int, real, real, ParabolicType, real);
void HlldMHD(DataBlock &, int, real, real);
void RoeMHD(DataBlock &, int, real, real);

/********************************
 * Local Kokkos Inline function *
 * ******************************/

KOKKOS_INLINE_FUNCTION void K_Flux(real &F_RHO, ARG_EXPAND(real &F_MX1, real &F_MX2, real &F_MX3, real &F_BX1, real &F_BX2, real &F_BX3, real &F_ENG),
                                   real V_RHO, ARG_EXPAND(real V_VX1, real V_VX2, real V_VX3, real V_BX1, real V_BX2, real V_BX3, real V_PRS),
                                   real U_RHO, ARG_EXPAND(real U_MX1, real U_MX2, real U_MX3, real U_BX1, real U_BX2, real U_BX3, real U_ENG),
                                   real C2Iso, int dir) {
    real U_MXn;
    EXPAND (
            real V_VXn;
            real V_BXn; ,
            real V_VXt;
            real V_BXt;  ,
            real V_VXb;
            real V_BXb;  )
    EXPAND (
        if (dir == IDIR) {
            U_MXn = U_MX1;
            EXPAND (
            V_VXn = V_VX1;
            V_BXn = V_BX1;  ,
            V_VXt = V_VX2;
            V_BXt = V_BX2;  ,
            V_VXb = V_VX3;
            V_BXb = V_BX3;  )
        }
    ,
        if (dir == JDIR) {
            U_MXn = U_MX2;
            EXPAND (
            V_VXn = V_VX2;
            V_BXn = V_BX2;  ,
            V_VXt = V_VX1;
            V_BXt = V_BX1;  ,
            V_VXb = V_VX3;
            V_BXb = V_BX3;  )
        }
    ,
        if (dir == KDIR) {
            U_MXn = U_MX3;
            EXPAND (
            V_VXn = V_VX3;
            V_BXn = V_BX3;  ,
            V_VXt = V_VX1;
            V_BXt = V_BX1;  ,
            V_VXb = V_VX2;
            V_BXb = V_BX2;  )
        }
    )
    
    F_RHO = U_MXn;
    EXPAND( F_MX1 = U_MX1*V_VXn - V_BXn*V_BX1; ,
            F_MX2 = U_MX2*V_VXn - V_BXn*V_BX2; ,
            F_MX3 = U_MX3*V_VXn - V_BXn*V_BX3;)

    EXPAND (
            real F_BXn;  ,
            real F_BXt;  ,
            real F_BXb;  )
    EXPAND(F_BXn = ZERO_F;                             ,
           F_BXt = V_VXn*V_BXt - V_BXn*V_VXt;   ,
           F_BXb = V_VXn*V_BXb - V_BXn*V_VXb; )

    real Bmag2 = EXPAND(V_BX1*V_BX1 , + V_BX2*V_BX2, + V_BX3*V_BX3);

    #if HAVE_ENERGY
    real ptot  = V_PRS + HALF_F*Bmag2;
    #elif defined(ISOTHERMAL)
    real ptot  = C2Iso * V_RHO + HALF_F*Bmag2;
    #else
    #error "K_Flux not defined for this EOS!"
    #endif

    #if HAVE_ENERGY
    F_ENG   = (U_ENG + ptot)*V_VXn - V_BXn * (EXPAND(V_VX1*V_BX1 , + V_VX2*V_BX2, + V_VX3*V_BX3));
    #endif

    // Add back pressure in the flux (not included in original PLUTO implementation)
    real F_MXn   = ptot;
    EXPAND (
        if (dir == IDIR) {
            F_MX1 += F_MXn;
            EXPAND(F_BX1 = F_BXn;  , 
                   F_BX2 = F_BXt;  , 
                   F_BX3 = F_BXb;  )
        }
    ,
        if (dir == JDIR) {
            F_MX2 += F_MXn;
            EXPAND(F_BX2 = F_BXn;  , 
                   F_BX1 = F_BXt;  , 
                   F_BX3 = F_BXb;  )
        }
    ,
        if (dir == KDIR) {
            F_MX3 += F_MXn;
            EXPAND(F_BX3 = F_BXn;  , 
                   F_BX1 = F_BXt;  , 
                   F_BX2 = F_BXb;  )
        }
    )
} 

KOKKOS_INLINE_FUNCTION void K_ConsToPrim(real &Vc_RHO, ARG_EXPAND(real &Vc_VX1, real &Vc_VX2, real &Vc_VX3, real &Vc_BX1, real &Vc_BX2, real &Vc_BX3, real &Vc_PRS),
                                         real Uc_RHO, ARG_EXPAND(real Uc_MX1, real Uc_MX2, real Uc_MX3, real Uc_BX1, real Uc_BX2, real Uc_BX3, real Uc_ENG),
                                         real gamma_m1) {
    

    Vc_RHO = Uc_RHO;

    EXPAND( Vc_VX1 = Uc_MX1/Uc_RHO; ,
            Vc_VX2 = Uc_MX2/Uc_RHO; ,
            Vc_VX3 = Uc_MX3/Uc_RHO;)

    EXPAND( Vc_BX1 = Uc_BX1;  ,
            Vc_BX2 = Uc_BX2;  ,
            Vc_BX3 = Uc_BX3;  )
         

#if HAVE_ENERGY
    real kin,mag;
    kin = HALF_F / Uc_RHO * (EXPAND(    Uc_MX1*Uc_MX1 , 
                                    + Uc_MX2*Uc_MX2 ,
                                    + Uc_MX3*Uc_MX3 ));

    mag = HALF_F * (EXPAND(    Uc_BX1*Uc_BX1 , 
                                    + Uc_BX2*Uc_BX2 ,
                                    + Uc_BX3*Uc_BX3));


    Vc_PRS = gamma_m1 * (Uc_ENG - kin - mag);
    
    // Check pressure positivity
    if(Vc_PRS<= ZERO_F) {
        #ifdef SMALL_PRESSURE_TEMPERATURE
        Vc_PRS = SMALL_PRESSURE_TEMPERATURE*Vc_RHO;
        #else
        Vc_PRS = SMALL_PRESSURE_FIX;
        #endif
        Uc_ENG = Vc_PRS/gamma_m1+kin+mag;
    }
#endif  // Have_energy
}

KOKKOS_INLINE_FUNCTION void K_PrimToCons(real &Uc_RHO, ARG_EXPAND(real &Uc_MX1, real &Uc_MX2, real &Uc_MX3, real &Uc_BX1, real &Uc_BX2, real &Uc_BX3, real &Uc_ENG),
                                         real Vc_RHO, ARG_EXPAND(real Vc_VX1, real Vc_VX2, real Vc_VX3, real Vc_BX1, real Vc_BX2, real Vc_BX3, real Vc_PRS),
                                         real gamma_m1) {

    Uc_RHO = Vc_RHO;

    EXPAND( Uc_MX1 = Vc_VX1*Vc_RHO; ,
            Uc_MX2 = Vc_VX2*Vc_RHO; ,
            Uc_MX3 = Vc_VX3*Vc_RHO;)


    EXPAND( Uc_BX1 = Vc_BX1;  ,
            Uc_BX2 = Vc_BX2;  ,
            Uc_BX3 = Vc_BX3;  )

#if HAVE_ENERGY

    Uc_ENG = Vc_PRS / gamma_m1 
                + HALF_F * Vc_RHO * (EXPAND(  Vc_VX1*Vc_VX1 , 
                                         + Vc_VX2*Vc_VX2 ,
                                         + Vc_VX3*Vc_VX3 ))
                + HALF_F * (EXPAND(    Uc_BX1*Uc_BX1 , 
                                     + Uc_BX2*Uc_BX2 ,
                                     + Uc_BX3*Uc_BX3)); 
#endif  // Have_energy

}
