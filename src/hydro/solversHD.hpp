/********************************
 * Local Kokkos Inline function *
 * ******************************/
#pragma once
#include "../idefix.hpp"

void TvdlfHD(DataBlock &, int, real, real);
void HllHD(DataBlock &, int, real, real);
void HllcHD(DataBlock &, int, real, real);
void RoeHD(DataBlock &, int, real, real);

static KOKKOS_INLINE_FUNCTION void K_Flux(real &F_RHO, ARG_EXPAND(real &F_MX1, real &F_MX2, real &F_MX3, real &F_ENG),
                                          const real V_RHO, ARG_EXPAND(const real V_VX1, const real V_VX2, const real V_VX3, const real V_PRS),
                                          const real U_RHO, ARG_EXPAND(const real U_MX1, const real U_MX2, const real U_MX3, const real U_ENG),
                                          const real C2Iso, const int dir) {
    real F_MXn, V_VXn;

#if HAVE_ENERGY
    F_MXn     = V_PRS;

#elif EOS == ISOTHERMAL
    // Add back pressure in the flux
    F_MXn     = C2Iso * V_RHO;
#endif
    
    EXPAND(F_MX1 = ZERO_F; ,
           F_MX2 = ZERO_F; ,
           F_MX3 = ZERO_F;)
    
    EXPAND (
        if (dir == IDIR) {
            F_RHO = U_MX1;
            F_MX1 = F_MXn;
            V_VXn = V_VX1;
        }
    ,
        if (dir == JDIR) {
            F_RHO = U_MX2;
            F_MX2 = F_MXn;
            V_VXn = V_VX2;
        }
    ,
        if (dir == KDIR) {
            F_RHO = U_MX3;
            F_MX3 = F_MXn;
            V_VXn = V_VX3;
        }
    )

    EXPAND(F_MX1 += U_MX1*V_VXn; ,
           F_MX2 += U_MX2*V_VXn; ,
           F_MX3 += U_MX3*V_VXn;)

    
#if HAVE_ENERGY
    F_ENG     = (U_ENG + V_PRS)*V_VXn;
#endif
    
} 

static KOKKOS_INLINE_FUNCTION void K_ConsToPrim(real &Vc_RHO, ARG_EXPAND(real &Vc_VX1, real &Vc_VX2, real &Vc_VX3, real &Vc_PRS),
                                                const real Uc_RHO, ARG_EXPAND(const real Uc_MX1, const real Uc_MX2, const real Uc_MX3, const real Uc_ENG),
                                                const real gamma_m1) {
    

    Vc_RHO = Uc_RHO;

    EXPAND( Vc_VX1 = Uc_MX1/Uc_RHO; ,
            Vc_VX2 = Uc_MX2/Uc_RHO; ,
            Vc_VX3 = Uc_MX3/Uc_RHO;)     

#if HAVE_ENERGY
    real kin;
    kin = HALF_F / Uc_RHO * (EXPAND(    Uc_MX1*Uc_MX1 , 
                                    + Uc_MX2*Uc_MX2 ,
                                    + Uc_MX3*Uc_MX3 ));

    Vc_PRS = gamma_m1 * (Uc_ENG - kin);
#endif  // Have_energy
}

static KOKKOS_INLINE_FUNCTION void K_PrimToCons(real &Uc_RHO, ARG_EXPAND(real &Uc_MX1, real &Uc_MX2, real &Uc_MX3, real &Uc_ENG),
                                                const real Vc_RHO, ARG_EXPAND(const real Vc_VX1, const real Vc_VX2, const real Vc_VX3, const real Vc_PRS),
                                                const real gamma_m1) {

    Uc_RHO = Vc_RHO;

    EXPAND( Uc_MX1 = Vc_VX1*Vc_RHO; ,
            Uc_MX2 = Vc_VX2*Vc_RHO; ,
            Uc_MX3 = Vc_VX3*Vc_RHO;)

#if HAVE_ENERGY

    Uc_ENG = Vc_PRS / gamma_m1 
                + HALF_F * Vc_RHO * EXPAND(  Vc_VX1*Vc_VX1 , 
                                         + Vc_VX2*Vc_VX2 ,
                                         + Vc_VX3*Vc_VX3 ); 
#endif  // Have_energy

}
