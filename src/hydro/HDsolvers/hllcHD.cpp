#include "../idefix.hpp"
#include "../solversHD.hpp"

// Compute Riemann fluxes from states using HLLC solver
void HllcHD(DataBlock & data, int dir, real gamma, real C2Iso) {
    int ioffset,joffset,koffset;

    idfx::pushRegion("HLLC_Solver");
    
    ioffset=joffset=koffset=0;
    // Determine the offset along which we do the extrapolation
    switch(dir) {
        case(IDIR): ioffset = 1;
            break;
        case(JDIR): joffset = 1;
            break;
        case(KDIR): koffset = 1;
            break;
        default:
            IDEFIX_ERROR("Wrong direction");
    }

    IdefixArray4D<real> PrimL = data.PrimL;
    IdefixArray4D<real> PrimR = data.PrimR;
    IdefixArray4D<real> Flux = data.FluxRiemann;
    IdefixArray1D<real> dx = data.dx[dir];
    IdefixArray3D<real> invDt = data.InvDtHyp;

    real gamma_m1 = gamma - ONE_F;

    idefix_for("HLLC_Kernel",data.beg[KDIR],data.end[KDIR]+koffset,data.beg[JDIR],data.end[JDIR]+joffset,data.beg[IDIR],data.end[IDIR]+ioffset,
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
                    
                    real vs;

                    EXPAND(
                    real uL_MXn; real uR_MXn;
                    real fluxL_MXn; real fluxR_MXn;
                    real usL_MXn; real usR_MXn;
                    ,
                    real vL_VXt; real vR_VXt;
                    real usL_MXt; real usR_MXt;
                    ,
                    real vL_VXb; real vR_VXb;
                    real uL_MXb; real uR_MXb;
                    );
                    
                    real usL_RHO;
                    EXPAND(
                    real usL_MX1; ,
                    real usL_MX2; ,
                    real usL_MX3; );
                
                    real usR_RHO;
                    EXPAND(
                    real usR_MX1; ,
                    real usR_MX2; ,
                    real usR_MX3; );
                    
                    EXPAND (
                    if (dir == IDIR) {
                        EXPAND (
                        uL_MXn = uL_MX1;
                        uR_MXn = uR_MX1;
                        fluxL_MXn = fluxL_MX1;
                        fluxR_MXn = fluxR_MX1;
                        ,
                        vL_VXt =  vL_VX2;
                        vR_VXt =  vR_VX2;
                        ,
                        vL_VXb =  vL_VX3;
                        vR_VXb =  vR_VX3;
                        )
                    }
                ,
                    if (dir == JDIR) {
                        EXPAND (
                        uL_MXn = uL_MX2;
                        uR_MXn = uR_MX2;
                        fluxL_MXn = fluxL_MX2;
                        fluxR_MXn = fluxR_MX2;
                        ,
                        vL_VXt =  vL_VX1;
                        vR_VXt =  vR_VX1;
                        ,
                        vL_VXb =  vL_VX3;
                        vR_VXb =  vR_VX3;
                        )
                    }
                ,
                    if (dir == KDIR) {
                        EXPAND (
                        uL_MXn = uL_MX3;
                        uR_MXn = uR_MX3;
                        fluxL_MXn = fluxL_MX3;
                        fluxR_MXn = fluxR_MX3;
                        ,
                        vL_VXt =  vL_VX1;
                        vR_VXt =  vR_VX1;
                        ,
                        vL_VXb =  vL_VX2;
                        vR_VXb =  vR_VX2;
                        )
                    }
                )
                
#if HAVE_ENERGY
                    real usL_ENG, usR_ENG;
                    real qL, qR, wL, wR;
                    qL = vL_PRS + uL_MXn*(vL_VXn - SL);
                    qR = vR_PRS + uR_MXn*(vR_VXn - SR);

                    wL = vL_RHO*(vL_VXn - SL);
                    wR = vR_RHO*(vR_VXn - SR);

                    vs = (qR - qL)/(wR - wL); // wR - wL > 0 since SL < 0, SR > 0

                    usL_RHO = uL_RHO*(SL - vL_VXn)/(SL - vs);
                    usR_RHO = uR_RHO*(SR - vR_VXn)/(SR - vs);
                    EXPAND(usL_MXn = usL_RHO*vs;     usR_MXn = usR_RHO*vs;      ,
                           usL_MXt = usL_RHO*vL_VXt; usR_MXt = usR_RHO*vR_VXt;  ,
                           usL_MXb = usL_RHO*vL_VXb; usR_MXb = usR_RHO*vR_VXb;)

                    usL_ENG =    uL_ENG/vL_RHO
                               + (vs - vL_VXn)*(vs + vL_PRS/(vL_RHO*(SL - vL_VXn)));
                    usR_ENG =    uR_ENG/vR_RHO
                               + (vs - vR_VXn)*(vs + vR_PRS/(vR_RHO*(SR - vR_VXn)));

                    usL_ENG *= usL_RHO;
                    usR_ENG *= usR_RHO;
#else
                    real scrh = 1.0/(SR - SL);
                    real rho  = (SR*uR_RHO - SL*uL_RHO - fluxR_RHO + fluxL_RHO)*scrh;
                    real mx   = (SR*uR_MXn - SL*uL_MXn - fluxR_MXn + fluxL_MXn)*scrh;

                    vs  = (  SR*fluxL_RHO - SL*fluxR_RHO
                           + SR*SL*(uR_RHO - uL_RHO));
                    vs *= scrh;
                    vs /= rho;
                    
                    usL_RHO = usR_RHO = rho;
                    EXPAND(usL_MXn = usR_MXn = mx;                     ,
                           usL_MXt = rho*vL_VXt; usR_MXt = rho*vR_VXt; ,
                           usL_MXb = rho*vL_VXb; usR_MXb = rho*vR_VXb;)
#endif
                    
                    EXPAND (
                        if (dir == IDIR) {
                            EXPAND (
                            usL_MX1 = usL_MXn;
                            usR_MX1 = usR_MXn;
                            ,
                            usL_MX2 = usL_MXt;
                            usR_MX2 = usR_MXt;
                            ,
                            usL_MX3 = usL_MXb;
                            usR_MX3 = usR_MXb;
                            )
                        }
                    ,
                        if (dir == JDIR) {
                            EXPAND (
                            usL_MX1 = usL_MXt;
                            usR_MX1 = usR_MXt;
                            ,
                            usL_MX2 = usL_MXn;
                            usR_MX2 = usR_MXn;
                            ,
                            usL_MX3 = usL_MXb;
                            usR_MX3 = usR_MXb;
                            )
                        }
                    ,
                        if (dir == KDIR) {
                            EXPAND (
                            usL_MX1 = usL_MXt;
                            usR_MX1 = usR_MXt;
                            ,
                            usL_MX2 = usL_MXb;
                            usR_MX2 = usR_MXb;
                            ,
                            usL_MX3 = usL_MXn;
                            usR_MX3 = usR_MXn;
                            )
                        }
                    )
                    
                    // 5-- Compute the flux from the left and right states
                    if (vs >= ZERO_F) {
                        Flux(RHO,k,j,i) = fluxL_RHO + SL*(usL_RHO - uL_RHO);
                        EXPAND (
                        Flux(MX1,k,j,i) = fluxL_MX1 + SL*(usL_MX1 - uL_MX1); ,
                        Flux(MX2,k,j,i) = fluxL_MX2 + SL*(usL_MX2 - uL_MX2); ,
                        Flux(MX3,k,j,i) = fluxL_MX3 + SL*(usL_MX3 - uL_MX3);
                        )
#if HAVE_ENERGY
                        Flux(ENG,k,j,i) = fluxL_ENG + SL*(usL_ENG - uL_ENG);
#endif
                    } else {
                        Flux(RHO,k,j,i) = fluxR_RHO + SR*(usR_RHO - uR_RHO);
                        EXPAND (
                        Flux(MX1,k,j,i) = fluxR_MX1 + SR*(usR_MX1 - uR_MX1); ,
                        Flux(MX2,k,j,i) = fluxR_MX2 + SR*(usR_MX2 - uR_MX2); ,
                        Flux(MX3,k,j,i) = fluxR_MX3 + SR*(usR_MX3 - uR_MX3);
                        )
#if HAVE_ENERGY
                        Flux(ENG,k,j,i) = fluxR_ENG + SR*(usR_ENG - uR_ENG);
#endif
                    }
                }

                //6-- Compute maximum dt for this sweep
                const int ig = ioffset*i + joffset*j + koffset*k;

                invDt(k,j,i) += cmax/dx(ig);
            });

    idfx::popRegion();

}
