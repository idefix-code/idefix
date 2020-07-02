#include "../idefix.hpp"
#include "../solversHD.hpp"

// Compute Riemann fluxes from states using TVDLF solver
void TvdlfHD(DataBlock & data, int dir, real gamma, real C2Iso) {
    int ioffset,joffset,koffset;

    idfx::pushRegion("TVDLF_Solver");
    ioffset=joffset=koffset=0;
    // Determine the offset along which we do the extrapolation
    if(dir==IDIR) ioffset=1;
    if(dir==JDIR) joffset=1;
    if(dir==KDIR) koffset=1;

    IdefixArray4D<real> PrimL = data.PrimL;
    IdefixArray4D<real> PrimR = data.PrimR;
    IdefixArray4D<real> Flux = data.FluxRiemann;
    IdefixArray1D<real> dx = data.dx[dir];
    IdefixArray3D<real> invDt = data.InvDtHyp;

    real gamma_m1=gamma-ONE_F;

    idefix_for("TVDLF_Kernel",data.beg[KDIR],data.end[KDIR]+koffset,data.beg[JDIR],data.end[JDIR]+joffset,data.beg[IDIR],data.end[IDIR]+ioffset,
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
                real cRL, cmax;

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
                
                // 1.3-- compute primitive variables in the center
                real vRL_VXn;
                EXPAND (
                    if (dir == IDIR) {
                        vRL_VXn = HALF_F * (vL_VX1 + vR_VX1);
                    }
                ,
                    if (dir == JDIR) {
                        vRL_VXn = HALF_F * (vL_VX2 + vR_VX2);
                    }
                ,
                    if (dir == KDIR) {
                        vRL_VXn = HALF_F * (vL_VX3 + vR_VX3);
                    }
                )
#if HAVE_ENERGY
                real vRL_RHO = HALF_F * (vL_RHO + vR_RHO);
                real vL_PRS = PrimL(PRS,k,j,i);
                real vR_PRS = PrimR(PRS,k,j,i);
                real vRL_PRS = HALF_F * (vL_PRS + vR_PRS);
                
                real uL_ENG, uR_ENG;
                real fluxL_ENG, fluxR_ENG;
#endif
                

                // 2-- Get the wave speed
#if HAVE_ENERGY
                cRL = SQRT((gamma_m1+ONE_F)*(vRL_PRS/vRL_RHO));
#else
                cRL = SQRT(C2Iso);
#endif
                cmax = FMAX(FABS(vRL_VXn+cRL),FABS(vRL_VXn-cRL));
                
                
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
                Flux(RHO,k,j,i) = HALF_F*(fluxL_RHO+fluxR_RHO - cmax*(uR_RHO-uL_RHO));
                EXPAND ( Flux(MX1,k,j,i) = HALF_F * (fluxL_MX1 + fluxR_MX1 - cmax*(uR_MX1-uL_MX1)); ,
                         Flux(MX2,k,j,i) = HALF_F * (fluxL_MX2 + fluxR_MX2 - cmax*(uR_MX2-uL_MX2)); ,
                         Flux(MX3,k,j,i) = HALF_F * (fluxL_MX3 + fluxR_MX3 - cmax*(uR_MX3-uL_MX3));
                       )
#if HAVE_ENERGY
                Flux(ENG,k,j,i) = HALF_F*(fluxL_ENG+fluxR_ENG - cmax*(uR_ENG-uL_ENG));
#endif
                
                //6-- Compute maximum dt for this sweep
                const int ig = ioffset*i + joffset*j + koffset*k;

                invDt(k,j,i) += cmax/dx(ig);

            });

    idfx::popRegion();

}

