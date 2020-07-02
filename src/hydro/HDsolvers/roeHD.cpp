#include "../idefix.hpp"
#include "../solversHD.hpp"

#define ROE_AVERAGE 0

// Compute Riemann fluxes from states using ROE solver
void RoeHD(DataBlock & data, int dir, real gamma, real C2Iso) {
    int ioffset,joffset,koffset;

    idfx::pushRegion("ROE_Solver");
    ioffset=joffset=koffset=0;
    // Determine the offset along which we do the extrapolation
    switch(dir) {
        case(IDIR):
            ioffset = 1;
            break;
        case(JDIR):
            joffset=1;
            break;
        case(KDIR):
            koffset=1;
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
#if HAVE_ENERGY
    real gmm1_inv = ONE_F / gamma_m1;
#endif
    
    idefix_for("ROE_Kernel",data.beg[KDIR],data.end[KDIR]+koffset,data.beg[JDIR],data.end[JDIR]+joffset,data.beg[IDIR],data.end[IDIR]+ioffset,
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
                
                real dv_RHO;
                EXPAND(
                    real dv_VX1; ,
                    real dv_VX2; ,
                    real dv_VX3; );

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

                // 1.2-- Compute dV
                dv_RHO = vR_RHO - vL_RHO;
                EXPAND(
                    dv_VX1 = vR_VX1 - vL_VX1; ,
                    dv_VX2 = vR_VX2 - vL_VX2; ,
                    dv_VX3 = vR_VX3 - vL_VX3; );
                
#if HAVE_ENERGY
                real vL_PRS = PrimL(PRS,k,j,i);
                real vR_PRS = PrimR(PRS,k,j,i);
                real dv_PRS = vR_PRS - vL_PRS;
                real uL_ENG, uR_ENG;
                real fluxL_ENG, fluxR_ENG;
#endif
                
                // 2-- Compute the conservative variables
                K_PrimToCons(uL_RHO, ARG_EXPAND(uL_MX1, uL_MX2, uL_MX3, uL_ENG),
                             vL_RHO, ARG_EXPAND(vL_VX1, vL_VX2, vL_VX3, vL_PRS),
                             gamma_m1);
                
                K_PrimToCons(uR_RHO, ARG_EXPAND(uR_MX1, uR_MX2, uR_MX3, uR_ENG),
                             vR_RHO, ARG_EXPAND(vR_VX1, vR_VX2, vR_VX3, vR_PRS),
                             gamma_m1);

                
                // 3-- Compute the left and right fluxes
                K_Flux(fluxL_RHO, ARG_EXPAND(fluxL_MX1, fluxL_MX2, fluxL_MX3, fluxL_ENG),
                       vL_RHO, ARG_EXPAND(vL_VX1, vL_VX2, vL_VX3, vL_PRS),
                       uL_RHO, ARG_EXPAND(uL_MX1, uL_MX2, uL_MX3, uL_ENG),
                       C2Iso, dir);
                
                K_Flux(fluxR_RHO, ARG_EXPAND(fluxR_MX1, fluxR_MX2, fluxR_MX3, fluxR_ENG),
                       vR_RHO, ARG_EXPAND(vR_VX1, vR_VX2, vR_VX3, vR_PRS),
                       uR_RHO, ARG_EXPAND(uR_MX1, uR_MX2, uR_MX3, uR_ENG),
                       C2Iso, dir);
                
                // 4-- Compute the square of the sound speed
                real a, a2, a2L, a2R;
#if HAVE_ENERGY
                a2L = gamma * vL_PRS / vL_RHO;
                a2R = gamma * vR_PRS / vR_RHO;
                
#else
                a2L = C2Iso;
                a2R = C2Iso;
#endif
                
                //  ----  Define Wave Jumps  ----
                real um_RHO;
                EXPAND(
                real um_VX1; ,
                real um_VX2; ,
                real um_VX3; );
#if HAVE_ENERGY
                real um_PRS;
                real h, hl, hr,vel2;
#endif
                
#if ROE_AVERAGE == YES
                real s, c;
                s       = sqrt(vR_RHO/vL_RHO);
                um_RHO  = vL_RHO * s;
                s       = ONE_F/(ONE_F + s); 
                c       = ONE_F - s;

                EXPAND(um_VX1 = s*vL_VX1 + c*vR_VX1;  ,
                    um_VX2 = s*vL_VX2 + c*vR_VX2;  ,
                    um_VX3 = s*vL_VX3 + c*vR_VX3;)

    #if HAVE_ENERGY
                
                vel2 = EXPAND(um_VX1*um_VX1, + um_VX2*um_VX2, + um_VX3*um_VX3);

                hl  = HALF_F*(EXPAND(vL_VX1*vL_VX1, + vL_VX2*vL_VX2, + vL_VX3*vL_VX3));    
                hl += a2L*gmm1_inv;

                hr = HALF_F*(EXPAND(vR_VX1*vR_VX1, + vR_VX2*vR_VX2, + vR_VX3*vR_VX3));    
                hr += a2R*gmm1_inv;

                h = s*hl + c*hr;

                /* -------------------------------------------------
                the following should be  equivalent to 

                scrh = EXPAND(   dv_VX1*dv_VX1,
                + dv_VX2*dv_VX2,
                + dv_VX3*dv_VX3);

                a2 = s*a2L + c*a2R + HALF_F*gamma_m1*s*c*scrh;

                and therefore always positive.
                just work out the coefficiendnts...
                -------------------------------------------------- */

                a2 = gamma_m1*(h - HALF_F*vel2);
                a  = sqrt(a2);
    #else
                a2 = HALF_F*(a2L + a2R);
                a  = sqrt(a2);
    #endif // HAVE_ENERGY
#else
                um_RHO = HALF_F*(vR_RHO+vL_RHO);
                EXPAND(
                    um_VX1 = HALF_F*(vR_VX1+vL_VX1); ,
                    um_VX2 = HALF_F*(vR_VX2+vL_VX2); ,
                    um_VX3 = HALF_F*(vR_VX3+vL_VX3); );

    #if HAVE_ENERGY
                um_PRS = HALF_F*(vR_PRS+vL_PRS);
                
                a2   = gamma*um_PRS/um_RHO;
                a    = sqrt(a2);

                vel2 = EXPAND(um_VX1*um_VX1, + um_VX2*um_VX2, + um_VX3*um_VX3);
                h    = HALF_F*vel2 + a2/gamma_m1;
    #else
                a2 = HALF_F*(a2L + a2R);
                a  = sqrt(a2);
    #endif // HAVE_ENERGY
#endif // ROE_AVERAGE

// **********************************************************************************
                /* ----------------------------------------------------------------
                define non-zero components of conservative eigenvectors Rc, 
                eigenvalues (lambda) and wave strenght eta = L.du     
                ----------------------------------------------------------------  */
                
                real lambda_RHO;
                EXPAND(
                    real lambda_MX1; ,
                    real lambda_MX2; ,
                    real lambda_MX3; )
                
                real eta_RHO;
                EXPAND(
                    real eta_MX1; ,
                    real eta_MX2; ,
                    real eta_MX3; )
                
#if HAVE_ENERGY
                real lambda_ENG;
                real eta_ENG;
#endif
                
                real Rc_RHO_RHO = ZERO_F;
                EXPAND(
                    real Rc_RHO_MX1 = ZERO_F; ,
                    real Rc_RHO_MX2 = ZERO_F; ,
                    real Rc_RHO_MX3 = ZERO_F; );
                
                EXPAND(
                    real Rc_MX1_RHO = ZERO_F;
                    EXPAND(
                        real Rc_MX1_MX1 = ZERO_F; ,
                        real Rc_MX1_MX2 = ZERO_F; ,
                        real Rc_MX1_MX3 = ZERO_F; );
                    ,
                    real Rc_MX2_RHO = ZERO_F;
                    EXPAND(
                        real Rc_MX2_MX1 = ZERO_F; ,
                        real Rc_MX2_MX2 = ZERO_F; ,
                        real Rc_MX2_MX3 = ZERO_F; );
                    ,
                    real Rc_MX3_RHO = ZERO_F;
                    EXPAND(
                        real Rc_MX3_MX1 = ZERO_F; ,
                        real Rc_MX3_MX2 = ZERO_F; ,
                        real Rc_MX3_MX3 = ZERO_F; );
                );
#if HAVE_ENERGY
                real Rc_RHO_ENG = ZERO_F;
                real Rc_MX1_ENG = ZERO_F;
                real Rc_MX2_ENG = ZERO_F;
                real Rc_MX3_ENG = ZERO_F;
                
                real Rc_ENG_RHO = ZERO_F;
                EXPAND(
                    real Rc_ENG_MX1 = ZERO_F; ,
                    real Rc_ENG_MX2 = ZERO_F; ,
                    real Rc_ENG_MX3 = ZERO_F; );
                real Rc_ENG_ENG = ZERO_F;
#endif
                
                EXPAND(
                    real um_VXn;
                    real dv_VXn;
                    real vL_VXn;
                    real vR_VXn;
                    real Rc_MXn_RHO;
                    real Rc_MXn_MX1;
                    ,
                    real um_VXt;
                    real dv_VXt;
                    real Rc_MXt_RHO;
                    real Rc_MXt_MX1;
                    real Rc_MXt_MX2;
                    ,
                    real um_VXb;
                    real dv_VXb;
                    real Rc_MXb_RHO;
                    real Rc_MXb_MX1;
                )
                
                EXPAND (
                    if (dir == IDIR) {
                        EXPAND (
                        um_VXn = um_VX1;
                        dv_VXn = dv_VX1;
                        vL_VXn = vL_VX1;
                        vR_VXn = vR_VX1;
                        ,
                        um_VXt =  um_VX2;
                        dv_VXt =  dv_VX2;
                        ,
                        um_VXb =  um_VX3;
                        dv_VXb =  dv_VX3;
                        )
                    }
                ,
                    if (dir == JDIR) {
                        EXPAND (
                        um_VXn = um_VX2;
                        dv_VXn = dv_VX2;
                        vL_VXn = vL_VX2;
                        vR_VXn = vR_VX2;
                        ,
                        um_VXt =  um_VX1;
                        dv_VXt =  dv_VX1;
                        ,
                        um_VXb =  um_VX3;
                        dv_VXb =  dv_VX3;
                        )
                    }
                ,
                    if (dir == KDIR) {
                        EXPAND (
                        um_VXn = um_VX3;
                        dv_VXn = dv_VX3;
                        vL_VXn = vL_VX3;
                        vR_VXn = vR_VX3;
                        ,
                        um_VXt =  um_VX1;
                        dv_VXt =  dv_VX1;
                        ,
                        um_VXb =  um_VX2;
                        dv_VXb =  dv_VX2;
                        )
                    }
                )
                
                //  ---- (u - c_s)  ---- 

                lambda_RHO = um_VXn - a;
#if HAVE_ENERGY
                eta_RHO = HALF_F/a2*(dv_PRS - dv_VXn*um_RHO*a);
#else
                eta_RHO = HALF_F*(dv_RHO - um_RHO*dv_VXn/a);
#endif

                Rc_RHO_RHO        = ONE_F;
                EXPAND(Rc_MXn_RHO = um_VXn - a;   ,
                Rc_MXt_RHO = um_VXt;              ,
                Rc_MXb_RHO = um_VXb;)
#if HAVE_ENERGY
                Rc_ENG_RHO = h - um_VXn*a;
#endif

                /*  ---- (u + c_s)  ----  */ 

                lambda_MX1 = um_VXn + a;
#if HAVE_ENERGY
                eta_MX1    = HALF_F/a2*(dv_PRS + dv_VXn*um_RHO*a);
#else
                eta_MX1 = HALF_F*(dv_RHO + um_RHO*dv_VXn/a);
#endif

                Rc_RHO_MX1        = ONE_F;
                EXPAND(Rc_MXn_MX1 = um_VXn + a;   ,
                Rc_MXt_MX1 = um_VXt;              ,
                Rc_MXb_MX1 = um_VXb;)
#if HAVE_ENERGY
                Rc_ENG_MX1 = h + um_VXn*a;
#endif

#if HAVE_ENERGY
                /*  ----  (u)  ----  */ 

                //ENG         = 2;
                lambda_ENG = um_VXn;
                eta_ENG    = dv_RHO - dv_PRS/a2;
                Rc_RHO_ENG        = ONE_F;
                EXPAND(Rc_MX1_ENG = um_VX1;   ,
                Rc_MX2_ENG = um_VX2;          ,
                Rc_MX3_ENG = um_VX3;)
                Rc_ENG_ENG        = HALF_F*vel2;
#endif

#if COMPONENTS > 1
                /*  ----  (u)  ----  */ 

                lambda_MX2 = um_VXn;
                eta_MX2    = um_RHO*dv_VXt;
                Rc_MXt_MX2 = ONE_F;
    #if HAVE_ENERGY
                Rc_ENG_MX2 = um_VXt;  
    #endif
#endif

#if COMPONENTS > 2
                /*  ----  (u)  ----  */ 

                lambda_MX3 = um_VXn;
                eta_MX3    = um_RHO*dv_VXb;
                Rc_MXb_MX3 = ONE_F;
    #if HAVE_ENERGY
                Rc_ENG_MX3 = um_VXb;  
    #endif
#endif

                real cmax = FABS(um_VXn) + a;

                /* ---------------------------------------------
                use the HLL flux function if the interface 
                lies within a strong shock.
                The effect of this switch is visible
                in the Mach reflection test.
                --------------------------------------------- */

                real scrh, scrh1;
                real bmin, bmax;
#if HAVE_ENERGY
                scrh  = FABS(vL_PRS - vR_PRS);
                scrh /= FMIN(vL_PRS,vR_PRS);
#else
                scrh  = FABS(vL_RHO - vR_RHO);
                scrh /= FMIN(vL_RHO,vR_RHO);
                scrh *= a*a;
#endif

#if DIMENSIONS > 1                
                if (scrh > HALF_F && (vR_VXn < vL_VXn)) {   /* -- tunable parameter -- */

                    bmin = FMIN(ZERO_F, lambda_RHO);
                    bmax = FMAX(ZERO_F, lambda_MX1);
                    scrh1 = ONE_F/(bmax - bmin);
                    
                    Flux(RHO,k,j,i) = (bmin*bmax*(uR_RHO - uL_RHO) + bmax*fluxL_RHO - bmin*fluxR_RHO)*scrh1;
                    EXPAND (
                    Flux(MX1,k,j,i) = (bmin*bmax*(uR_MX1 - uL_MX1) + bmax*fluxL_MX1 - bmin*fluxR_MX1)*scrh1; ,
                    Flux(MX2,k,j,i) = (bmin*bmax*(uR_MX2 - uL_MX2) + bmax*fluxL_MX2 - bmin*fluxR_MX2)*scrh1; ,
                    Flux(MX3,k,j,i) = (bmin*bmax*(uR_MX3 - uL_MX3) + bmax*fluxL_MX3 - bmin*fluxR_MX3)*scrh1;
                    )
    #if HAVE_ENERGY
                Flux(ENG,k,j,i) = (bmin*bmax*(uR_ENG - uL_ENG) + bmax*fluxL_ENG - bmin*fluxR_ENG)*scrh1;
    #endif
                }
                else {
#endif
                    /* -----------------------------------------------------------
                                        compute Roe flux 
                    ----------------------------------------------------------- */

                    EXPAND (
                        if (dir == IDIR) {
                            EXPAND (
                            Rc_MX1_RHO = Rc_MXn_RHO;
                            Rc_MX1_MX1 = Rc_MXn_MX1;
                            ,
                            Rc_MX2_RHO = Rc_MXt_RHO;
                            Rc_MX2_MX1 = Rc_MXt_MX1;
                            Rc_MX2_MX2 = Rc_MXt_MX2;
                            ,
                            Rc_MX3_RHO = Rc_MXb_RHO;
                            Rc_MX3_MX1 = Rc_MXb_MX1;
                            Rc_MX3_MX3 = Rc_MXb_MX3;
                            )
                        }
                    ,
                        if (dir == JDIR) {
                            EXPAND (
                            Rc_MX1_RHO = Rc_MXt_RHO;
                            Rc_MX1_MX1 = Rc_MXt_MX1;
                            Rc_MX1_MX2 = Rc_MXt_MX2;
                            ,
                            Rc_MX2_RHO = Rc_MXn_RHO;
                            Rc_MX2_MX1 = Rc_MXn_MX1;
                            ,
                            Rc_MX3_RHO = Rc_MXb_RHO;
                            Rc_MX3_MX1 = Rc_MXb_MX1;
                            Rc_MX3_MX3 = Rc_MXb_MX3;
                            )
                        }
                    ,
                        if (dir == KDIR) {
                            EXPAND (
                            Rc_MX1_RHO = Rc_MXt_RHO;
                            Rc_MX1_MX1 = Rc_MXt_MX1;
                            Rc_MX1_MX2 = Rc_MXt_MX2;
                            ,
                            Rc_MX2_RHO = Rc_MXb_RHO;
                            Rc_MX2_MX1 = Rc_MXb_MX1;
                            Rc_MX2_MX3 = Rc_MXb_MX3;
                            ,
                            Rc_MX3_RHO = Rc_MXn_RHO;
                            Rc_MX3_MX1 = Rc_MXn_MX1;
                            )
                        }
                    )
                    
                    real alambda_RHO = FABS(lambda_RHO);
                    EXPAND(
                    real alambda_MX1 = FABS(lambda_MX1); ,
                    real alambda_MX2 = FABS(lambda_MX2); ,
                    real alambda_MX3 = FABS(lambda_MX3); )
#if HAVE_ENERGY
                    real alambda_ENG = FABS(lambda_ENG);
#endif

                    /*  ----  entropy fix  ----  */
                    real delta = 1.e-7;
                    if (alambda_RHO <= delta) {
                        alambda_RHO = HALF_F*lambda_RHO*lambda_RHO/delta + HALF_F*delta;
                    }
                    if (alambda_MX1 <= delta) {
                        alambda_MX1 = HALF_F*lambda_MX1*lambda_MX1/delta + HALF_F*delta;
                    }

                    alambda_RHO *= eta_RHO;
                    EXPAND(
                    alambda_MX1 *= eta_MX1; ,
                    alambda_MX2 *= eta_MX2; ,
                    alambda_MX3 *= eta_MX3; )
#if HAVE_ENERGY
                    alambda_ENG *= eta_ENG;
#endif
                    
                    real fl_RHO = alambda_RHO * Rc_RHO_RHO 
                        EXPAND(+ alambda_MX1 * Rc_RHO_MX1,
                               + alambda_MX2 * Rc_RHO_MX2,
                               + alambda_MX3 * Rc_RHO_MX3);
                    EXPAND(
                    real fl_MX1 = alambda_RHO * Rc_MX1_RHO
                        EXPAND(+ alambda_MX1 * Rc_MX1_MX1,
                               + alambda_MX2 * Rc_MX1_MX2,
                               + alambda_MX3 * Rc_MX1_MX3); ,
                    real fl_MX2 = alambda_RHO * Rc_MX2_RHO
                        EXPAND(+ alambda_MX1 * Rc_MX2_MX1,
                               + alambda_MX2 * Rc_MX2_MX2,
                               + alambda_MX3 * Rc_MX2_MX3); ,
                    real fl_MX3 = alambda_RHO * Rc_MX3_RHO
                        EXPAND(+ alambda_MX1 * Rc_MX3_MX1,
                               + alambda_MX2 * Rc_MX3_MX2,
                               + alambda_MX3 * Rc_MX3_MX3); )
#if HAVE_ENERGY
                    real fl_ENG = alambda_RHO * Rc_ENG_RHO
                        EXPAND(+ alambda_MX1 * Rc_ENG_MX1,
                               + alambda_MX2 * Rc_ENG_MX2,
                               + alambda_MX3 * Rc_ENG_MX3)
                               + alambda_ENG * Rc_ENG_ENG;
                    
                    fl_RHO += alambda_ENG * Rc_RHO_ENG;
                    EXPAND(
                    fl_MX1 += alambda_ENG * Rc_MX1_ENG; ,
                    fl_MX2 += alambda_ENG * Rc_MX2_ENG; ,
                    fl_MX3 += alambda_ENG * Rc_MX3_ENG; )
#endif
                    
                    Flux(RHO,k,j,i) = (fluxL_RHO + fluxR_RHO - fl_RHO)*HALF_F;
                    EXPAND (
                    Flux(MX1,k,j,i) = (fluxL_MX1 + fluxR_MX1 - fl_MX1)*HALF_F; ,
                    Flux(MX2,k,j,i) = (fluxL_MX2 + fluxR_MX2 - fl_MX2)*HALF_F; ,
                    Flux(MX3,k,j,i) = (fluxL_MX3 + fluxR_MX3 - fl_MX3)*HALF_F;
                    )
#if HAVE_ENERGY
                    Flux(ENG,k,j,i) = (fluxL_ENG + fluxR_ENG - fl_ENG)*HALF_F;
#endif
                    
#if DIMENSIONS > 1
                }
#endif

                //6-- Compute maximum dt for this sweep
                const int ig = ioffset*i + joffset*j + koffset*k;

                invDt(k,j,i) += cmax/dx(ig);

            });

    idfx::popRegion();

}
