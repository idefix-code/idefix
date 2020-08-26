#include "../idefix.hpp"
#include "solversHD.hpp"

// Compute Riemann fluxes from states using HLL solver
void HllHD(DataBlock & data, int dir, real gamma, real C2Iso) {
    int ioffset,joffset,koffset;

    idfx::pushRegion("HLL_Solver");
    ioffset=joffset=koffset=0;
    // Determine the offset along which we do the extrapolation
    if(dir==IDIR) ioffset=1;
    if(dir==JDIR) joffset=1;
    if(dir==KDIR) koffset=1;

    IdefixArray4D<real> PrimL = data.PrimL;
    IdefixArray4D<real> PrimR = data.PrimR;
    IdefixArray4D<real> Flux = data.FluxRiemann;
    IdefixArray1D<real> dx = data.dx[dir];
    IdefixArray3D<real> cMax = data.cMax;

    real gamma_m1 = gamma - ONE_F;

    idefix_for("HLL_Kernel",data.beg[KDIR],data.end[KDIR]+koffset,data.beg[JDIR],data.end[JDIR]+joffset,data.beg[IDIR],data.end[IDIR]+ioffset,
                        KOKKOS_LAMBDA (int k, int j, int i) 
            {

                // Primitive variables
                real vL_RHO;
                EXPAND(
                    real vL_VX1; ,
                    real vL_VX2; ,
                    real vL_VX3; );
                
                real vR_RHO;
                EXPAND(
                    real vR_VX1; ,
                    real vR_VX2; ,
                    real vR_VX3; );

                // Conservative variables
                real uL_RHO;
                EXPAND(
                    real uL_MX1; ,
                    real uL_MX2; ,
                    real uL_MX3; );
                
                real uR_RHO;
                EXPAND(
                    real uR_MX1; ,
                    real uR_MX2; ,
                    real uR_MX3; );

                // Flux (left and right)
                real fluxL_RHO;
                EXPAND(
                    real fluxL_MX1; ,
                    real fluxL_MX2; ,
                    real fluxL_MX3; );
                
                real fluxR_RHO;
                EXPAND(
                    real fluxR_MX1; ,
                    real fluxR_MX2; ,
                    real fluxR_MX3; );

                // Signal speeds
                real cL, cR, cmax;

                
                // 1.1-- Read primitive variables on the left state
                vL_RHO = PrimL(RHO,k,j,i);
                EXPAND(
                    vL_VX1 = PrimL(MX1,k,j,i); ,
                    vL_VX2 = PrimL(MX2,k,j,i); ,
                    vL_VX3 = PrimL(MX3,k,j,i); );
                
                // 1.2-- Read primitive variables on the right state
                vR_RHO = PrimR(RHO,k,j,i);
                EXPAND(
                    vR_VX1 = PrimR(MX1,k,j,i); ,
                    vR_VX2 = PrimR(MX2,k,j,i); ,
                    vR_VX3 = PrimR(MX3,k,j,i); );

                
                // 2-- Get the wave speed
#if HAVE_ENERGY
                real uL_ENG, uR_ENG;
                real fluxL_ENG, fluxR_ENG;
                
                real vL_PRS = PrimL(PRS,k,j,i);
                real vR_PRS = PrimR(PRS,k,j,i);
                
                cL = SQRT( gamma *(vL_PRS/vL_RHO));
                cR = SQRT( gamma *(vR_PRS/vR_RHO));
#else
                cL = SQRT(C2Iso);
                cR = cL;
#endif
                
                real vL_VXn, vR_VXn;
                
                EXPAND (
                    if (dir == IDIR) {
                        vL_VXn = vL_VX1;
                        vR_VXn = vR_VX1;
                    }
                ,
                    if (dir == JDIR) {
                        vL_VXn = vL_VX2;
                        vR_VXn = vR_VX2;
                    }
                ,
                    if (dir == KDIR) {
                        vL_VXn = vL_VX3;
                        vR_VXn = vR_VX3;
                    }
                )
                
                real cminL = vL_VXn - cL;
                real cmaxL = vL_VXn + cL;
                
                real cminR = vR_VXn - cR;
                real cmaxR = vR_VXn + cR;
                
                real SL = FMIN(cminL, cminR);
                real SR = FMAX(cmaxL, cmaxR);
                
                cmax  = FMAX(FABS(SL), FABS(SR));
                
                
                // 3-- Compute the conservative variables
                K_PrimToCons(uL_RHO, ARG_EXPAND(uL_MX1, uL_MX2, uL_MX3, uL_ENG),
                             vL_RHO, ARG_EXPAND(vL_VX1, vL_VX2, vL_VX3, vL_PRS),
                             gamma_m1);
                
                K_PrimToCons(uR_RHO, ARG_EXPAND(uR_MX1, uR_MX2, uR_MX3, uR_ENG),
                             vR_RHO, ARG_EXPAND(vR_VX1, vR_VX2, vR_VX3, vR_PRS),
                             gamma_m1);

                
                // 4-- Compute the left and right fluxes
                K_Flux(fluxL_RHO, ARG_EXPAND(fluxL_MX1, fluxL_MX2, fluxL_MX3, fluxL_ENG),
                       vL_RHO, ARG_EXPAND(vL_VX1, vL_VX2, vL_VX3, vL_PRS),
                       uL_RHO, ARG_EXPAND(uL_MX1, uL_MX2, uL_MX3, uL_ENG),
                       C2Iso, dir);
                
                K_Flux(fluxR_RHO, ARG_EXPAND(fluxR_MX1, fluxR_MX2, fluxR_MX3, fluxR_ENG),
                       vR_RHO, ARG_EXPAND(vR_VX1, vR_VX2, vR_VX3, vR_PRS),
                       uR_RHO, ARG_EXPAND(uR_MX1, uR_MX2, uR_MX3, uR_ENG),
                       C2Iso, dir);

                
                // 5-- Compute the flux from the left and right states
                if (SL > 0) {
                    
                    Flux(RHO,k,j,i) = fluxL_RHO;
                    EXPAND ( Flux(MX1,k,j,i) = fluxL_MX1; ,
                            Flux(MX2,k,j,i) = fluxL_MX2; ,
                            Flux(MX3,k,j,i) = fluxL_MX3;
                        )
#if HAVE_ENERGY
                    Flux(ENG,k,j,i) = fluxL_ENG;
#endif
                    
                }
                else if (SR < 0) {
                    
                    Flux(RHO,k,j,i) = fluxR_RHO;
                    EXPAND ( Flux(MX1,k,j,i) = fluxR_MX1; ,
                            Flux(MX2,k,j,i) = fluxR_MX2; ,
                            Flux(MX3,k,j,i) = fluxR_MX3;
                        )
#if HAVE_ENERGY
                    Flux(ENG,k,j,i) = fluxR_ENG;
#endif
                    
                }
                else {
                    
                    Flux(RHO,k,j,i) = (SL*SR*(uR_RHO - uL_RHO) + SR*fluxL_RHO - SL*fluxR_RHO) / (SR - SL) ;
                    EXPAND ( Flux(MX1,k,j,i) = (SL*SR*(uR_MX1 - uL_MX1) + SR*fluxL_MX1 - SL*fluxR_MX1) / (SR - SL); ,
                            Flux(MX2,k,j,i) = (SL*SR*(uR_MX2 - uL_MX2) + SR*fluxL_MX2 - SL*fluxR_MX2) / (SR - SL); ,
                            Flux(MX3,k,j,i) = (SL*SR*(uR_MX3 - uL_MX3) + SR*fluxL_MX3 - SL*fluxR_MX3) / (SR - SL);
                        )
#if HAVE_ENERGY
                    Flux(ENG,k,j,i) = (SL*SR*(uR_ENG - uL_ENG) + SR*fluxL_ENG - SL*fluxR_ENG) / (SR - SL);
#endif
                    
                }


                //6-- Compute maximum wave speed for this sweep

                cMax(k,j,i) = cmax;

            });

    idfx::popRegion();

}
