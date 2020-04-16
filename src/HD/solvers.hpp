/********************************
 * Local Kokkos Inline function *
 * ******************************/
#pragma once
#include "../idefix.hpp"

void Tvdlf(DataBlock &, int, real, real);
void Hll(DataBlock &, int, real, real);
void Hllc(DataBlock &, int, real, real);
void Roe(DataBlock &, int, real, real);

static KOKKOS_INLINE_FUNCTION void K_Flux(real *KOKKOS_RESTRICT F, const real *KOKKOS_RESTRICT V, const real *KOKKOS_RESTRICT U, real C2Iso, const int dir) {
    int VXn = VX1+dir;
    int MXn = VXn;

    F[RHO] = U[VXn];

    EXPAND(F[MX1] = U[MX1]*V[VXn]; ,
           F[MX2] = U[MX2]*V[VXn]; ,
           F[MX3] = U[MX3]*V[VXn];)


#if HAVE_ENERGY
    F[ENG]     = (U[ENG] + V[PRS])*V[VXn];
    F[MXn]     += V[PRS];

#elif EOS == ISOTHERMAL
    // Add back pressure in the flux
    F[MXn]     += C2Iso * V[RHO];
#endif
} 

static KOKKOS_INLINE_FUNCTION void K_ConsToPrim(real *KOKKOS_RESTRICT Vc, const real *KOKKOS_RESTRICT Uc, real gamma_m1) {
    

    Vc[RHO] = Uc[RHO];

    EXPAND( Vc[VX1] = Uc[MX1]/Uc[RHO]; ,
            Vc[VX2] = Uc[MX2]/Uc[RHO]; ,
            Vc[VX3] = Uc[MX3]/Uc[RHO];)     

#if HAVE_ENERGY
    real kin;
    kin = 0.5 / Uc[RHO] * (EXPAND(    Uc[MX1]*Uc[MX1] , 
                                    + Uc[MX2]*Uc[MX2] ,
                                    + Uc[MX3]*Uc[MX3] ));

    Vc[PRS] = gamma_m1 * (Uc[ENG] - kin);
#endif  // Have_energy
}

static KOKKOS_INLINE_FUNCTION void K_PrimToCons(real *KOKKOS_RESTRICT Uc, const real *KOKKOS_RESTRICT Vc, real gamma_m1) {

    Uc[RHO] = Vc[RHO];

    EXPAND( Uc[MX1] = Vc[VX1]*Vc[RHO]; ,
            Uc[MX2] = Vc[VX2]*Vc[RHO]; ,
            Uc[MX3] = Vc[VX3]*Vc[RHO];)

#if HAVE_ENERGY

    Uc[ENG] = Vc[PRS] / gamma_m1 
                + 0.5 * Vc[RHO] * EXPAND(  Vc[VX1]*Vc[VX1] , 
                                         + Vc[VX2]*Vc[VX2] ,
                                         + Vc[VX3]*Vc[VX3] ); 
#endif  // Have_energy

}
