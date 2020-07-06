#include "../idefix.hpp"
#include "solversMHD.hpp"

// Compute Riemann fluxes from states using HLLD solver
void HlldMHD(DataBlock & data, int dir, real gamma, real C2Iso) {

    int ioffset,joffset,koffset;
    int iextend, jextend,kextend;

    idfx::pushRegion("HLLD_MHD");
    
    ioffset=joffset=koffset=0;
    // extension in perp to the direction of integration, as required by CT.
    iextend=jextend=kextend=0;

    IdefixArray4D<real> PrimL = data.PrimL;
    IdefixArray4D<real> PrimR = data.PrimR;
    IdefixArray4D<real> Flux = data.FluxRiemann;
    IdefixArray1D<real> dx = data.dx[dir];
    IdefixArray3D<real> invDt = data.InvDtHyp;

    // References to required emf components
    IdefixArray3D<real> Eb;
    IdefixArray3D<real> Et;
    
    IdefixArray3D<int> SV;


    real gamma_m1=gamma-ONE_F;

    // Define normal, tangent and bi-tanget indices

    EXPAND(       ,
        int BXt;  ,
        int BXb;  )
#if !HAVE_ENERGY
    EXPAND(
        int BXn; int MXn;  ,
        int MXt;           ,
        int MXb;           )
#endif

    // st and sb will be useful only when Hall is included
    D_EXPAND( real st;  ,
                        ,
              real sb;  )
    
    switch(dir) {
        case(IDIR):
            ioffset = 1;
            D_EXPAND(
                    st = -ONE_F;  ,
                    jextend = 1;  ,
                    kextend = 1;
                    sb = +ONE_F;  )
#if !HAVE_ENERGY
            EXPAND( MXn = MX1; 
                    BXn = BX1;    , 
                    MXt = MX2;    , 
                    MXb = MX3;    )
#endif
            EXPAND(               , 
                    BXt = BX2;    ,
                    BXb = BX3;    )

            Et = data.emf.ezi;
            Eb = data.emf.eyi;
            SV = data.emf.svx;

            break;
        case(JDIR):
            joffset=1;
            D_EXPAND(
                    iextend = 1;
                    st = +ONE_F;  ,
                                  ,
                    kextend = 1;
                    sb = -ONE_F;  )
#if !HAVE_ENERGY
            EXPAND(MXn = MX2; 
                   BXn = BX2;     , 
                   MXt = MX1;     , 
                   MXb = MX3;     )
#endif
            EXPAND(               , 
                   BXt = BX1;     , 
                   BXb = BX3;     )

            Et = data.emf.ezj;
            Eb = data.emf.exj;
            SV = data.emf.svy;
            
            break;
        case(KDIR):
            koffset=1;
            D_EXPAND(
                    iextend = 1;
                    st = -ONE_F;  ,
                    jextend = 1;  ,
                    sb = +ONE_F;  )
#if !HAVE_ENERGY
            EXPAND(MXn = MX3; 
                   BXn = BX3;     , 
                   MXt = MX1;     , 
                   MXb = MX2;     )
#endif
            EXPAND(               , 
                   BXt = BX1;     , 
                   BXb = BX2;     )

            Et = data.emf.eyk;
            Eb = data.emf.exk;
            SV = data.emf.svz;
            
            break;
        default:
            IDEFIX_ERROR("Wrong direction");
    }

    idefix_for("CalcRiemannFlux",data.beg[KDIR]-kextend,data.end[KDIR]+koffset+kextend,data.beg[JDIR]-jextend,data.end[JDIR]+joffset+jextend,data.beg[IDIR]-iextend,data.end[IDIR]+ioffset+iextend,
                KOKKOS_LAMBDA (int k, int j, int i) 
    {
        // Primitive variables
        real vL_RHO;
        EXPAND(
            real vL_VX1;
            real vL_BX1; ,
            real vL_VX2;
            real vL_BX2; ,
            real vL_VX3;
            real vL_BX3; )
        
        real vR_RHO;
        EXPAND(
            real vR_VX1;
            real vR_BX1; ,
            real vR_VX2;
            real vR_BX2; ,
            real vR_VX3;
            real vR_BX3; )

        // Conservative variables
        real uL_RHO;
        EXPAND(
            real uL_MX1; 
            real uL_BX1; ,
            real uL_MX2;
            real uL_BX2; ,
            real uL_MX3;
            real uL_BX3; )
        
        real uR_RHO;
        EXPAND(
            real uR_MX1; 
            real uR_BX1; ,
            real uR_MX2;
            real uR_BX2; ,
            real uR_MX3;
            real uR_BX3; )
        
        // Flux (left and right)
        real fluxL_RHO;
        EXPAND(
            real fluxL_MX1; 
            real fluxL_BX1; ,
            real fluxL_MX2;
            real fluxL_BX2; ,
            real fluxL_MX3;
            real fluxL_BX3; )
        
        real fluxR_RHO;
        EXPAND(
            real fluxR_MX1; 
            real fluxR_BX1; ,
            real fluxR_MX2;
            real fluxR_BX2; ,
            real fluxR_MX3;
            real fluxR_BX3; )

        // Signal speeds
        real cL, cR, cmax;

        // 1.1-- Read primitive variables on the left state
        vL_RHO = PrimL(RHO,k,j,i);
        EXPAND(
            vL_VX1 = PrimL(VX1,k,j,i);
            vL_BX1 = PrimL(BX1,k,j,i); ,
            vL_VX2 = PrimL(VX2,k,j,i);
            vL_BX2 = PrimL(BX2,k,j,i); ,
            vL_VX3 = PrimL(VX3,k,j,i);
            vL_BX3 = PrimL(BX3,k,j,i); )
        
        // 1.2-- Read primitive variables on the right state
        vR_RHO = PrimR(RHO,k,j,i);
        EXPAND(
            vR_VX1 = PrimR(VX1,k,j,i);
            vR_BX1 = PrimR(BX1,k,j,i); ,
            vR_VX2 = PrimR(VX2,k,j,i);
            vR_BX2 = PrimR(BX2,k,j,i); ,
            vR_VX3 = PrimR(VX3,k,j,i);
            vR_BX3 = PrimR(BX3,k,j,i); )
        
#if HAVE_ENERGY
        real vL_PRS = PrimL(PRS,k,j,i);
        real vR_PRS = PrimR(PRS,k,j,i);
        real uL_ENG, uR_ENG;
        real fluxL_ENG, fluxR_ENG;
#endif

        EXPAND(
            real vL_VXn; real vR_VXn;
            real vL_BXn; real vR_BXn; ,
            real vL_BXt; real vR_BXt; ,
            real vL_BXb; real vR_BXb; )
        
        EXPAND (
        if (dir == IDIR) {
            EXPAND (
            vL_VXn = vL_VX1;
            vR_VXn = vR_VX1;
            vL_BXn = vL_BX1;
            vR_BXn = vR_BX1;  ,
            vL_BXt = vL_BX2;
            vR_BXt = vR_BX2;  ,
            vL_BXb = vL_BX3;
            vR_BXb = vR_BX3;  )
        }
        ,
        if (dir == JDIR) {
            EXPAND (
            vL_VXn = vL_VX2;
            vR_VXn = vR_VX2;
            vL_BXn = vL_BX2;
            vR_BXn = vR_BX2;  ,
            vL_BXt = vL_BX1;
            vR_BXt = vR_BX1;  ,
            vL_BXb = vL_BX3;
            vR_BXb = vR_BX3;  )
        }
        ,
        if (dir == KDIR) {
            EXPAND (
            vL_VXn = vL_VX3;
            vR_VXn = vR_VX3;
            vL_BXn = vL_BX3;
            vR_BXn = vR_BX3;  ,
            vL_BXt = vL_BX1;
            vR_BXt = vR_BX1;  ,
            vL_BXb = vL_BX2;
            vR_BXb = vR_BX2;  )
        }
        )
            
        // 2-- Get the wave speed
        real gpr, b1, b2, b3, Btmag2, Bmag2;
#if HAVE_ENERGY
        gpr = gamma*vL_PRS;
#else
        gpr = C2Iso*vL_RHO;
#endif

        // -- get total field
        b1 = b2 = b3 = ZERO_F;
        EXPAND (b1 = vL_BXn;  ,
                b2 = vL_BXt;  ,
                b3 = vL_BXb;)

        Btmag2 = b2*b2 + b3*b3;
        Bmag2  = b1*b1 + Btmag2;

        cL = gpr - Bmag2;
        cL = gpr + Bmag2 + sqrt(cL*cL + 4.0*gpr*Btmag2);
        cL = sqrt(HALF_F*cL/vL_RHO);
        
#if HAVE_ENERGY
        gpr = gamma*vR_PRS;
#else
        gpr = C2Iso*vR_RHO;
#endif

        // -- get total field
        b1 = b2 = b3 = ZERO_F;
        EXPAND (b1 = vR_BXn;  ,
                b2 = vR_BXt;  ,
                b3 = vR_BXb;)

        Btmag2 = b2*b2 + b3*b3;
        Bmag2  = b1*b1 + Btmag2;

        cR = gpr - Bmag2;
        cR = gpr + Bmag2 + sqrt(cR*cR + 4.0*gpr*Btmag2);
        cR = sqrt(HALF_F*cR/vR_RHO);
        
        // 4.1 
        real cminL = vL_VXn - cL;
        real cmaxL = vL_VXn + cL;
        
        real cminR = vR_VXn - cR;
        real cmaxR = vR_VXn + cR;
        
        real SL = FMIN(cminL, cminR);
        real SR = FMAX(cmaxL, cmaxR);
        
        cmax  = FMAX(FABS(SL), FABS(SR));

        // 2-- Compute the conservative variables
        K_PrimToCons(uL_RHO, ARG_EXPAND(uL_MX1, uL_MX2, uL_MX3, uL_BX1, uL_BX2, uL_BX3, uL_ENG),
                     vL_RHO, ARG_EXPAND(vL_VX1, vL_VX2, vL_VX3, vL_BX1, vL_BX2, vL_BX3, vL_PRS),
                     gamma_m1);
        K_PrimToCons(uR_RHO, ARG_EXPAND(uR_MX1, uR_MX2, uR_MX3, uR_BX1, uR_BX2, uR_BX3, uR_ENG),
                     vR_RHO, ARG_EXPAND(vR_VX1, vR_VX2, vR_VX3, vR_BX1, vR_BX2, vR_BX3, vR_PRS),
                     gamma_m1);

        // 3-- Compute the left and right fluxes
        K_Flux(fluxL_RHO, ARG_EXPAND(fluxL_MX1, fluxL_MX2, fluxL_MX3, fluxL_BX1, fluxL_BX2, fluxL_BX3, fluxL_ENG),
               vL_RHO, ARG_EXPAND(vL_VX1, vL_VX2, vL_VX3, vL_BX1, vL_BX2, vL_BX3, vL_PRS),
               uL_RHO, ARG_EXPAND(uL_MX1, uL_MX2, uL_MX3, uL_BX1, uL_BX2, uL_BX3, uL_ENG),
               C2Iso, dir);
        K_Flux(fluxR_RHO, ARG_EXPAND(fluxR_MX1, fluxR_MX2, fluxR_MX3, fluxR_BX1, fluxR_BX2, fluxR_BX3, fluxR_ENG),
               vR_RHO, ARG_EXPAND(vR_VX1, vR_VX2, vR_VX3, vR_BX1, vR_BX2, vR_BX3, vR_PRS),
               uR_RHO, ARG_EXPAND(uR_MX1, uR_MX2, uR_MX3, uR_BX1, uR_BX2, uR_BX3, uR_ENG),
               C2Iso, dir);

        // 5-- Compute the flux from the left and right states
        if (SL > 0){
            Flux(RHO,k,j,i) = fluxL_RHO;
            EXPAND ( Flux(MX1,k,j,i) = fluxL_MX1;
                 Flux(BX1,k,j,i) = fluxL_BX1;    ,
                 Flux(MX2,k,j,i) = fluxL_MX2;
                 Flux(BX2,k,j,i) = fluxL_BX2;    ,
                 Flux(MX3,k,j,i) = fluxL_MX3;
                 Flux(BX3,k,j,i) = fluxL_BX3;    )
#if HAVE_ENERGY
            Flux(ENG,k,j,i) = fluxL_ENG;
#endif
        }
        else if (SR < 0) {
            Flux(RHO,k,j,i) = fluxR_RHO;
            EXPAND ( Flux(MX1,k,j,i) = fluxR_MX1;
                 Flux(BX1,k,j,i) = fluxR_BX1;    ,
                 Flux(MX2,k,j,i) = fluxR_MX2;
                 Flux(BX2,k,j,i) = fluxR_BX2;    ,
                 Flux(MX3,k,j,i) = fluxR_MX3;
                 Flux(BX3,k,j,i) = fluxR_BX3;    )
#if HAVE_ENERGY
            Flux(ENG,k,j,i) = fluxR_ENG;
#endif
        }
        else {
            
            EXPAND(                   ,
                real vL_VXt; real vR_VXt; ,
                real vL_VXb; real vR_VXb; )
            
            EXPAND (
            if (dir == IDIR) {
                EXPAND (   ,
                vL_VXt = vL_VX2;
                vR_VXt = vR_VX2;  ,
                vL_VXb = vL_VX3;
                vR_VXb = vR_VX3;  )
            }
            ,
            if (dir == JDIR) {
                EXPAND (   ,
                vL_VXt = vL_VX1;
                vR_VXt = vR_VX1;  ,
                vL_VXb = vL_VX3;
                vR_VXb = vR_VX3;  )
            }
            ,
            if (dir == KDIR) {
                EXPAND (   ,
                vL_VXt = vL_VX1;
                vR_VXt = vR_VX1;  ,
                vL_VXb = vL_VX2;
                vR_VXb = vR_VX2;  )
            }
            )
            
            EXPAND(
                real uL_MXn; real uR_MXn;  ,
                real uL_BXt; real uR_BXt;  ,
                real uL_BXb; real uR_BXb;  )
            
            EXPAND (
            if (dir == IDIR) {
                EXPAND (
                uL_MXn = uL_MX1;
                uR_MXn = uR_MX1;  ,
                uL_BXt = uL_BX2;
                uR_BXt = uR_BX2;  ,
                uL_BXb = uL_BX3;
                uR_BXb = uR_BX3;  )
            }
            ,
            if (dir == JDIR) {
                EXPAND (
                uL_MXn = uL_MX2;
                uR_MXn = uR_MX2;  ,
                uL_BXt = uL_BX1;
                uR_BXt = uR_BX1;  ,
                uL_BXb = uL_BX3;
                uR_BXb = uR_BX3;  )
            }
            ,
            if (dir == KDIR) {
                EXPAND (
                uL_MXn = uL_MX3;
                uR_MXn = uR_MX3;  ,
                uL_BXt = uL_BX1;
                uR_BXt = uR_BX1;  ,
                uL_BXb = uL_BX2;
                uR_BXb = uR_BX2;  )
            }
            )
            
            real scrh, scrhL, scrhR, duL, duR, sBx, Bx, SM, S1L, S1R;
            
#if HAVE_ENERGY
            
            real pts, sqrL, sqrR, Bx1;
            EXPAND(                  ,
                real vsL; real vsR;  ,
                real wsL; real wsR;  )
            int revert_to_hllc;
            
            real ptL  = vL_PRS + HALF_F* ( EXPAND(vL_BX1*vL_BX1 , + vL_BX2*vL_BX2, + vL_BX3*vL_BX3) );
            real ptR  = vR_PRS + HALF_F* ( EXPAND(vR_BX1*vR_BX1 , + vR_BX2*vR_BX2, + vR_BX3*vR_BX3) );
            
            real usL_RHO;
            EXPAND(
                real usL_MX1;
                real usL_BX1; ,
                real usL_MX2;
                real usL_BX2; ,
                real usL_MX3;
                real usL_BX3; )
        
            real usR_RHO;
            EXPAND(
                real usR_MX1;
                real usR_BX1; ,
                real usR_MX2;
                real usR_BX2; ,
                real usR_MX3;
                real usR_BX3; )

            real usL_ENG, usR_ENG;
            
            EXPAND(
                real usL_MXn; real usR_MXn;
                real usL_BXn; real usR_BXn;  ,
                real usL_MXt; real usR_MXt;
                real usL_BXt; real usR_BXt;  ,
                real usL_MXb; real usR_MXb;
                real usL_BXb; real usR_BXb;  )
            
            // 3c. Compute U*(L), U^*(R)

            scrh = ONE_F/(SR - SL);
            Bx1  = Bx = (SR*vR_BXn - SL*vL_BXn)*scrh; 
            sBx  = (Bx > ZERO_F ? ONE_F : -ONE_F);

            duL  = SL - vL_VXn;
            duR  = SR - vR_VXn;

            scrh = ONE_F/(duR*uR_RHO - duL*uL_RHO);
            SM   = (duR*uR_MXn - duL*uL_MXn - ptR + ptL)*scrh;

            pts  = duR*uR_RHO*ptL - duL*uL_RHO*ptR + 
                    vL_RHO*vR_RHO*duR*duL*(vR_VXn- vL_VXn);
            pts *= scrh;

            usL_RHO = uL_RHO*duL/(SL - SM);
            usR_RHO = uR_RHO*duR/(SR - SM);

            sqrL = sqrt(usL_RHO);
            sqrR = sqrt(usR_RHO);

            S1L = SM - fabs(Bx)/sqrL;
            S1R = SM + fabs(Bx)/sqrR;

            // -----------------------------------------------------------------
            // 3d When S1L -> SL or S1R -> SR a degeneracy occurs. 
            // Although Miyoshi & Kusano say that no jump exists, we don't
            // think this is actually true. 
            // Indeed, vy*, vz*, By*, Bz* cannot be solved independently. 
            // In this case we revert to the HLLC solver of Li (2005),  except
            // for the term v.B in the region, which we compute in our own way.
            // Note, that by comparing the expressions of Li (2005) and
            // Miyoshi & Kusano (2005), the only change involves a 
            // re-definition of By* and Bz* in terms of By(HLL), Bz(HLL).
            // -----------------------------------------------------------------

            revert_to_hllc = 0;

            if ( (S1L - SL) <  1.e-4*(SM - SL) ) revert_to_hllc = 1;
            if ( (S1R - SR) > -1.e-4*(SR - SM) ) revert_to_hllc = 1;

            if (revert_to_hllc) {

                EXPAND(
                    real Uhll_BX1; ,
                    real Uhll_BX2; ,
                    real Uhll_BX3; )
                
                //Uhll_RHO = (SR*uR_RHO - SL*uL_RHO + fluxL_RHO - fluxR_RHO) / (SR - SL);
                EXPAND ( Uhll_BX1 = (SR*uR_BX1 - SL*uL_BX1 + fluxL_BX1 - fluxR_BX1) / (SR - SL);   ,
                         Uhll_BX2 = (SR*uR_BX2 - SL*uL_BX2 + fluxL_BX2 - fluxR_BX2) / (SR - SL);   ,
                         Uhll_BX3 = (SR*uR_BX3 - SL*uL_BX3 + fluxL_BX3 - fluxR_BX3) / (SR - SL);   )
                
                EXPAND(
                    real Uhll_BXn;  ,
                    real Uhll_BXt;  ,
                    real Uhll_BXb;  )
            
                EXPAND (
                if (dir == IDIR) {
                    EXPAND (
                    Uhll_BXn = Uhll_BX1;  ,
                    Uhll_BXt = Uhll_BX2;  ,
                    Uhll_BXb = Uhll_BX3;  )
                }
                ,
                if (dir == JDIR) {
                    EXPAND (
                    Uhll_BXn = Uhll_BX2;  ,
                    Uhll_BXt = Uhll_BX1;  ,
                    Uhll_BXb = Uhll_BX3;  )
                }
                ,
                if (dir == KDIR) {
                    EXPAND (
                    Uhll_BXn = Uhll_BX3;  ,
                    Uhll_BXt = Uhll_BX1;  ,
                    Uhll_BXb = Uhll_BX2;  )
                }
                )
                
                EXPAND(usL_BXn = usR_BXn = Uhll_BXn;   ,
                        usL_BXt = usR_BXt = Uhll_BXt;  ,
                        usL_BXb = usR_BXb = Uhll_BXb;  )
                
                S1L = S1R = SM; // region ** should never be computed since
                                // fluxes are given in terms of UL* and UR*
            }
            else {
                // 3e. Compute states in the * regions

                scrhL = (uL_RHO*duL*duL - Bx*Bx)/(uL_RHO*duL*(SL - SM) - Bx*Bx);
                scrhR = (uR_RHO*duR*duR - Bx*Bx)/(uR_RHO*duR*(SR - SM) - Bx*Bx);
        
                EXPAND(usL_BXn  = Bx1;        ,
                    usL_BXt  = uL_BXt*scrhL;  ,
                    usL_BXb  = uL_BXb*scrhL;  )

                EXPAND(usR_BXn = Bx1;         ,
                    usR_BXt = uR_BXt*scrhR;   ,   
                    usR_BXb = uR_BXb*scrhR;   )

            }

            scrhL = Bx/(uL_RHO*duL);
            scrhR = Bx/(uR_RHO*duR);

            EXPAND(                                          ;  ,
                    vsL = vL_VXt - scrhL*(usL_BXt - uL_BXt);
                    vsR = vR_VXt - scrhR*(usR_BXt - uR_BXt);  ,

                    wsL = vL_VXb - scrhL*(usL_BXb - uL_BXb);
                    wsR = vR_VXb - scrhR*(usR_BXb - uR_BXb); )

            EXPAND(usL_MXn = usL_RHO*SM; 
                    usR_MXn = usR_RHO*SM;   ,
            
                    usL_MXt = usL_RHO*vsL;
                    usR_MXt = usR_RHO*vsR;  ,

                    usL_MXb = usL_RHO*wsL;
                    usR_MXb = usR_RHO*wsR;)

            // -- Energy --

            scrhL  = EXPAND(vL_VXn*Bx1, + vL_VXt*uL_BXt, + vL_VXb*uL_BXb);
            scrhL -= EXPAND(     SM*Bx1, +    vsL*usL_BXt, +    wsL*usL_BXb);
            usL_ENG  = duL*uL_ENG - ptL*vL_VXn + pts*SM + Bx*scrhL;
            usL_ENG /= SL - SM;

            scrhR  = EXPAND(vR_VXn*Bx1, + vR_VXt*uR_BXt, + vR_VXb*uR_BXb);
            scrhR -= EXPAND(     SM*Bx1, +    vsR*usR_BXt, +    wsR*usR_BXb);
            usR_ENG = duR*uR_ENG - ptR*vR_VXn + pts*SM + Bx*scrhR;
            usR_ENG /= SR - SM;

            EXPAND (
                if (dir == IDIR) {
                    EXPAND (
                    usL_MX1 = usL_MXn;
                    usL_BX1 = usL_BXn; ,
                    usL_MX2 = usL_MXt;
                    usL_BX2 = usL_BXt; ,
                    usL_MX3 = usL_MXb;
                    usL_BX3 = usL_BXb; )
                }
                ,
                if (dir == JDIR) {
                    EXPAND (
                    usL_MX2 = usL_MXn;
                    usL_BX2 = usL_BXn; ,
                    usL_MX1 = usL_MXt;
                    usL_BX1 = usL_BXt; ,
                    usL_MX3 = usL_MXb;
                    usL_BX3 = usL_BXb; )
                }
                ,
                if (dir == KDIR) {
                    EXPAND (
                    usL_MX3 = usL_MXn;
                    usL_BX3 = usL_BXn; ,
                    usL_MX1 = usL_MXt;
                    usL_BX1 = usL_BXt; ,
                    usL_MX2 = usL_MXb;
                    usL_BX2 = usL_BXb; )
                }
            )
            
            EXPAND (
                if (dir == IDIR) {
                    EXPAND (
                    usR_MX1 = usR_MXn;
                    usR_BX1 = usR_BXn; ,
                    usR_MX2 = usR_MXt;
                    usR_BX2 = usR_BXt; ,
                    usR_MX3 = usR_MXb;
                    usR_BX3 = usR_BXb; )
                }
                ,
                if (dir == JDIR) {
                    EXPAND (
                    usR_MX2 = usR_MXn;
                    usR_BX2 = usR_BXn; ,
                    usR_MX1 = usR_MXt;
                    usR_BX1 = usR_BXt; ,
                    usR_MX3 = usR_MXb;
                    usR_BX3 = usR_BXb; )
                }
                ,
                if (dir == KDIR) {
                    EXPAND (
                    usR_MX3 = usR_MXn;
                    usR_BX3 = usR_BXn; ,
                    usR_MX1 = usR_MXt;
                    usR_BX1 = usR_BXt; ,
                    usR_MX2 = usR_MXb;
                    usR_BX2 = usR_BXb; )
                }
            )
            
        // 3c. Compute flux when S1L > 0 or S1R < 0

            if (S1L >= ZERO_F) {       //  ----  Region L*
                Flux(RHO,k,j,i) = fluxL_RHO + SL*(usL_RHO - uL_RHO);
                EXPAND ( Flux(MX1,k,j,i) = fluxL_MX1 + SL*(usL_MX1 - uL_MX1);
                        Flux(BX1,k,j,i) = fluxL_BX1 + SL*(usL_BX1 - uL_BX1);    ,
                        Flux(MX2,k,j,i) = fluxL_MX2 + SL*(usL_MX2 - uL_MX2);
                        Flux(BX2,k,j,i) = fluxL_BX2 + SL*(usL_BX2 - uL_BX2);    ,
                        Flux(MX3,k,j,i) = fluxL_MX3 + SL*(usL_MX3 - uL_MX3);
                        Flux(BX3,k,j,i) = fluxL_BX3 + SL*(usL_BX3 - uL_BX3);    )
                
                Flux(ENG,k,j,i) = fluxL_ENG + SL*(usL_ENG - uL_ENG);
                
            }
            else if (S1R <= ZERO_F) {    //  ----  Region R*
                Flux(RHO,k,j,i) = fluxR_RHO + SR*(usR_RHO - uR_RHO);
                EXPAND ( Flux(MX1,k,j,i) = fluxR_MX1 + SR*(usR_MX1 - uR_MX1);
                        Flux(BX1,k,j,i) = fluxR_BX1 + SR*(usR_BX1 - uR_BX1);    ,
                        Flux(MX2,k,j,i) = fluxR_MX2 + SR*(usR_MX2 - uR_MX2);
                        Flux(BX2,k,j,i) = fluxR_BX2 + SR*(usR_BX2 - uR_BX2);    ,
                        Flux(MX3,k,j,i) = fluxR_MX3 + SR*(usR_MX3 - uR_MX3);
                        Flux(BX3,k,j,i) = fluxR_BX3 + SR*(usR_BX3 - uR_BX3);    )
                
                Flux(ENG,k,j,i) = fluxR_ENG + SR*(usR_ENG - uR_ENG);
                
            } 
            else {   // -- This state exists only if B_x != 0

                // Compute U**
                EXPAND(        ,
                    real vss;  ,
                    real wss;  )
                
                real ussl_RHO;
                EXPAND(
                    real ussl_MX1;
                    real ussl_BX1; ,
                    real ussl_MX2;
                    real ussl_BX2; ,
                    real ussl_MX3;
                    real ussl_BX3; )
                real ussl_ENG;
                
                
                real ussr_RHO;
                EXPAND(
                    real ussr_MX1;
                    real ussr_BX1; ,
                    real ussr_MX2;
                    real ussr_BX2; ,
                    real ussr_MX3;
                    real ussr_BX3; )
                real ussr_ENG;
                
                EXPAND(real ussl_MXn;
                    real ussr_MXn;  ,
                    real ussl_MXt;
                    real ussr_MXt;  ,
                    real ussl_MXb;
                    real ussr_MXb;  )
                
                EXPAND(real ussl_BXn;
                    real ussr_BXn;  ,
                    real ussl_BXt;      
                    real ussr_BXt;  ,
                    real ussl_BXb;        
                    real ussr_BXb;  )

                ussl_RHO = usL_RHO;
                ussr_RHO = usR_RHO;
                
                EXPAND(                           ,
            
                    vss  = sqrL*vsL + sqrR*vsR + (usR_BXt - usL_BXt)*sBx;       
                    vss /= sqrL + sqrR;        ,
                    
                    wss  = sqrL*wsL + sqrR*wsR + (usR_BXb - usL_BXb)*sBx;
                    wss /= sqrL + sqrR;)

                EXPAND(ussl_MXn = ussl_RHO*SM;
                    ussr_MXn = ussr_RHO*SM;    ,
            
                    ussl_MXt = ussl_RHO*vss;
                    ussr_MXt = ussr_RHO*vss;  ,
                
                    ussl_MXb = ussl_RHO*wss;
                    ussr_MXb = ussr_RHO*wss;)           
            
                EXPAND(ussl_BXn = ussr_BXn = Bx1;   ,

                    ussl_BXt  = sqrL*usR_BXt + sqrR*usL_BXt + sqrL*sqrR*(vsR - vsL)*sBx;
                    ussl_BXt /= sqrL + sqrR;        
                    ussr_BXt  = ussl_BXt;        ,
                
                    ussl_BXb  = sqrL*usR_BXb + sqrR*usL_BXb + sqrL*sqrR*(wsR - wsL)*sBx;
                    ussl_BXb /= sqrL + sqrR;        
                    ussr_BXb  = ussl_BXb;)
                
                // -- Energy jump

                scrhL  = EXPAND(SM*Bx1, +  vsL*usL_BXt, +  wsL*usL_BXb);
                scrhL -= EXPAND(SM*Bx1, +  vss*ussl_BXt, +  wss*ussl_BXb);

                scrhR  = EXPAND(SM*Bx1, +  vsR*usR_BXt, +  wsR*usR_BXb);
                scrhR -= EXPAND(SM*Bx1, +  vss*ussr_BXt, +  wss*ussr_BXb);

                ussl_ENG = usL_ENG - sqrL*scrhL*sBx;
                ussr_ENG = usR_ENG + sqrR*scrhR*sBx;

                EXPAND (
                    if (dir == IDIR) {
                        EXPAND (
                        ussl_MX1 = ussl_MXn;
                        ussl_BX1 = ussl_BXn; ,
                        ussl_MX2 = ussl_MXt;
                        ussl_BX2 = ussl_BXt; ,
                        ussl_MX3 = ussl_MXb;
                        ussl_BX3 = ussl_BXb; )
                    }
                    ,
                    if (dir == JDIR) {
                        EXPAND (
                        ussl_MX2 = ussl_MXn;
                        ussl_BX2 = ussl_BXn; ,
                        ussl_MX1 = ussl_MXt;
                        ussl_BX1 = ussl_BXt; ,
                        ussl_MX3 = ussl_MXb;
                        ussl_BX3 = ussl_BXb; )
                    }
                    ,
                    if (dir == KDIR) {
                        EXPAND (
                        ussl_MX3 = ussl_MXn;
                        ussl_BX3 = ussl_BXn; ,
                        ussl_MX1 = ussl_MXt;
                        ussl_BX1 = ussl_BXt; ,
                        ussl_MX2 = ussl_MXb;
                        ussl_BX2 = ussl_BXb; )
                    }
                )
                
                EXPAND (
                    if (dir == IDIR) {
                        EXPAND (
                        ussr_MX1 = ussr_MXn;
                        ussr_BX1 = ussr_BXn; ,
                        ussr_MX2 = ussr_MXt;
                        ussr_BX2 = ussr_BXt; ,
                        ussr_MX3 = ussr_MXb;
                        ussr_BX3 = ussr_BXb; )
                    }
                    ,
                    if (dir == JDIR) {
                        EXPAND (
                        ussr_MX2 = ussr_MXn;
                        ussr_BX2 = ussr_BXn; ,
                        ussr_MX1 = ussr_MXt;
                        ussr_BX1 = ussr_BXt; ,
                        ussr_MX3 = ussr_MXb;
                        ussr_BX3 = ussr_BXb; )
                    }
                    ,
                    if (dir == KDIR) {
                        EXPAND (
                        ussr_MX3 = ussr_MXn;
                        ussr_BX3 = ussr_BXn; ,
                        ussr_MX1 = ussr_MXt;
                        ussr_BX1 = ussr_BXt; ,
                        ussr_MX2 = ussr_MXb;
                        ussr_BX2 = ussr_BXb; )
                    }
                )

                if (SM >= ZERO_F) {           //  ----  Region L**
                    Flux(RHO,k,j,i) = fluxL_RHO + S1L*(ussl_RHO - usL_RHO) + SL*(usL_RHO - uL_RHO);
                    EXPAND ( Flux(MX1,k,j,i) = fluxL_MX1 + S1L*(ussl_MX1 - usL_MX1) + SL*(usL_MX1 - uL_MX1);
                            Flux(BX1,k,j,i) = fluxL_BX1 + S1L*(ussl_BX1 - usL_BX1) + SL*(usL_BX1 - uL_BX1);    ,
                            Flux(MX2,k,j,i) = fluxL_MX2 + S1L*(ussl_MX2 - usL_MX2) + SL*(usL_MX2 - uL_MX2);
                            Flux(BX2,k,j,i) = fluxL_BX2 + S1L*(ussl_BX2 - usL_BX2) + SL*(usL_BX2 - uL_BX2);    ,
                            Flux(MX3,k,j,i) = fluxL_MX3 + S1L*(ussl_MX3 - usL_MX3) + SL*(usL_MX3 - uL_MX3);
                            Flux(BX3,k,j,i) = fluxL_BX3 + S1L*(ussl_BX3 - usL_BX3) + SL*(usL_BX3 - uL_BX3);    )
                
                    Flux(ENG,k,j,i) = fluxL_ENG + S1L*(ussl_ENG - usL_ENG) + SL*(usL_ENG - uL_ENG);

                }else{                   //  ----  Region R**
                    Flux(RHO,k,j,i) = fluxR_RHO + S1R*(ussr_RHO - usR_RHO) + SR*(usR_RHO - uR_RHO);
                    EXPAND ( Flux(MX1,k,j,i) = fluxR_MX1 + S1R*(ussr_MX1 - usR_MX1) + SR*(usR_MX1 - uR_MX1);
                            Flux(BX1,k,j,i) = fluxR_BX1 + S1R*(ussr_BX1 - usR_BX1) + SR*(usR_BX1 - uR_BX1);    ,
                            Flux(MX2,k,j,i) = fluxR_MX2 + S1R*(ussr_MX2 - usR_MX2) + SR*(usR_MX2 - uR_MX2);
                            Flux(BX2,k,j,i) = fluxR_BX2 + S1R*(ussr_BX2 - usR_BX2) + SR*(usR_BX2 - uR_BX2);    ,
                            Flux(MX3,k,j,i) = fluxR_MX3 + S1R*(ussr_MX3 - usR_MX3) + SR*(usR_MX3 - uR_MX3);
                            Flux(BX3,k,j,i) = fluxR_BX3 + S1R*(ussr_BX3 - usR_BX3) + SR*(usR_BX3 - uR_BX3);    )
                
                    Flux(ENG,k,j,i) = fluxR_ENG + S1R*(ussr_ENG - usR_ENG) + SR*(usR_ENG - uR_ENG);
                    
                }
            }  // end if (S1L < 0 S1R > 0)

#else
            real uL_BXn; real uR_BXn;
            
            EXPAND (
            if (dir == IDIR) {
                EXPAND (
                uL_MXn = uL_MX1;
                uR_MXn = uR_MX1;
                uL_BXn = uL_BX1;
                uR_BXn = uR_BX1;  ,
                uL_BXt = uL_BX2;
                uR_BXt = uR_BX2;  ,
                uL_BXb = uL_BX3;
                uR_BXb = uR_BX3;  )
            }
            ,
            if (dir == JDIR) {
                EXPAND (
                uL_MXn = uL_MX2;
                uR_MXn = uR_MX2;
                uL_BXn = uL_BX2;
                uR_BXn = uR_BX2;  ,
                uL_BXt = uL_BX1;
                uR_BXt = uR_BX1;  ,
                uL_BXb = uL_BX3;
                uR_BXb = uR_BX3;  )
            }
            ,
            if (dir == KDIR) {
                EXPAND (
                uL_MXn = uL_MX3;
                uR_MXn = uR_MX3;
                uL_BXn = uL_BX3;
                uR_BXn = uR_BX3;  ,
                uL_BXt = uL_BX1;
                uR_BXt = uR_BX1;  ,
                uL_BXb = uL_BX2;
                uR_BXb = uR_BX2;  )
            }
            )

            real rho, sqrho;
            int revert_to_hll;
            
            scrh = ONE_F/(SR - SL);
            duL = SL - vL_VXn;
            duR = SR - vR_VXn;

            Bx = (SR*vR_BXn - SL*vL_BXn)*scrh; 

            rho                = (uR_RHO*duR - uL_RHO*duL)*scrh;
            real Flux_RHO = (SL*uR_RHO*duR - SR*uL_RHO*duL)*scrh;
            
            //  compute S*
            
            sqrho = sqrt(rho);

            SM  = Flux_RHO/rho;
            S1L = SM - fabs(Bx)/sqrho;
            S1R = SM + fabs(Bx)/sqrho;

            // ---------------------------------------------
            //   Prevent degeneracies when S1L -> SL or 
            //   S1R -> SR. Revert to HLL if necessary.
            // ---------------------------------------------
            
            revert_to_hll = 0;

            if ( (S1L - SL) <  1.e-4*(SR - SL) ) revert_to_hll = 1;
            if ( (S1R - SR) > -1.e-4*(SR - SL) ) revert_to_hll = 1;

            if (revert_to_hll) {

                Flux(RHO,k,j,i) = (SL*SR*(uR_RHO - uL_RHO) + SR*fluxL_RHO - SL*fluxR_RHO) * scrh;
                EXPAND ( Flux(MX1,k,j,i) = (SL*SR*(uR_MX1 - uL_MX1) + SR*fluxL_MX1 - SL*fluxR_MX1) * scrh;
                        Flux(BX1,k,j,i) = (SL*SR*(uR_BX1 - uL_BX1) + SR*fluxL_BX1 - SL*fluxR_BX1) * scrh;    ,
                        Flux(MX2,k,j,i) = (SL*SR*(uR_MX2 - uL_MX2) + SR*fluxL_MX2 - SL*fluxR_MX2) * scrh;
                        Flux(BX2,k,j,i) = (SL*SR*(uR_BX2 - uL_BX2) + SR*fluxL_BX2 - SL*fluxR_BX2) * scrh;    ,
                        Flux(MX3,k,j,i) = (SL*SR*(uR_MX3 - uL_MX3) + SR*fluxL_MX3 - SL*fluxR_MX3) * scrh;
                        Flux(BX3,k,j,i) = (SL*SR*(uR_BX3 - uL_BX3) + SR*fluxL_BX3 - SL*fluxR_BX3) * scrh;    )

            }
            else {
                
                EXPAND(   ,
                    real uL_MXt; real uR_MXt; ,
                    real uL_MXb; real uR_MXb; )
                
                EXPAND (
                if (dir == IDIR) {
                    EXPAND (   ,
                    uL_MXt = uL_MX2;
                    uR_MXt = uR_MX2;  ,
                    uL_MXb = uL_MX3;
                    uR_MXb = uR_MX3;  )
                }
                ,
                if (dir == JDIR) {
                    EXPAND (   ,
                    uL_MXt = uL_MX1;
                    uR_MXt = uR_MX1;  ,
                    uL_MXb = uL_MX3;
                    uR_MXb = uR_MX3;  )
                }
                ,
                if (dir == KDIR) {
                    EXPAND (   ,
                    uL_MXt = uL_MX1;
                    uR_MXt = uR_MX1;  ,
                    uL_MXb = uL_MX2;
                    uR_MXb = uR_MX2;  )
                }
                )
                
                EXPAND(
                    real fluxL_MXn; real fluxR_MXn; ,
                    real fluxL_MXt; real fluxR_MXt;
                    real fluxL_BXt; real fluxR_BXt; ,
                    real fluxL_MXb; real fluxR_MXb;
                    real fluxL_BXb; real fluxR_BXb; )
                
                EXPAND (
                if (dir == IDIR) {
                    EXPAND (
                    fluxL_MXn = fluxL_MX1;
                    fluxR_MXn = fluxR_MX1;  ,
                    fluxL_MXt = fluxL_MX2;
                    fluxR_MXt = fluxR_MX2;
                    fluxL_BXt = fluxL_BX2;
                    fluxR_BXt = fluxR_BX2;  ,
                    fluxL_MXb = fluxL_MX3;
                    fluxR_MXb = fluxR_MX3;
                    fluxL_BXb = fluxL_BX3;
                    fluxR_BXb = fluxR_BX3;  )
                }
                ,
                if (dir == JDIR) {
                    EXPAND (
                    fluxL_MXn = fluxL_MX2;
                    fluxR_MXn = fluxR_MX2;  ,
                    fluxL_MXt = fluxL_MX1;
                    fluxR_MXt = fluxR_MX1;
                    fluxL_BXt = fluxL_BX1;
                    fluxR_BXt = fluxR_BX1;  ,
                    fluxL_MXb = fluxL_MX3;
                    fluxR_MXb = fluxR_MX3;
                    fluxL_BXb = fluxL_BX3;
                    fluxR_BXb = fluxR_BX3;  )
                }
                ,
                if (dir == KDIR) {
                    EXPAND (
                    fluxL_MXn = fluxL_MX3;
                    fluxR_MXn = fluxR_MX3;  ,
                    fluxL_MXt = fluxL_MX1;
                    fluxR_MXt = fluxR_MX1;
                    fluxL_BXt = fluxL_BX1;
                    fluxR_BXt = fluxR_BX1;  ,
                    fluxL_MXb = fluxL_MX2;
                    fluxR_MXb = fluxR_MX2;
                    fluxL_BXb = fluxL_BX2;
                    fluxR_BXb = fluxR_BX2;  )
                }
                )
                
                EXPAND(    ,
                    real usL_MXt; real usR_MXt;
                    real usL_BXt; real usR_BXt;  ,
                    real usL_MXb; real usR_MXb;
                    real usL_BXb; real usR_BXb;  )
                
                Flux(RHO,k,j,i) = Flux_RHO;
                Flux(MXn,k,j,i) = (SR*fluxL_MXn - SL*fluxR_MXn 
                                        + SR*SL*(uR_MXn - uL_MXn))*scrh;

                Flux(BXn,k,j,i) = SR*SL*(uR_BXn - uL_BXn)*scrh;


            //  Compute U*
                
                scrhL = ONE_F/((SL - S1L)*(SL - S1R));
                scrhR = ONE_F/((SR - S1L)*(SR - S1R));

                EXPAND(                                                        ,
                        usL_MXt = rho*vL_VXt - Bx*uL_BXt*(SM - vL_VXn)*scrhL;
                        usR_MXt = rho*vR_VXt - Bx*uR_BXt*(SM - vR_VXn)*scrhR;  ,
                        usL_MXb = rho*vL_VXb - Bx*uL_BXb*(SM - vL_VXn)*scrhL;
                        usR_MXb = rho*vR_VXb - Bx*uR_BXb*(SM - vR_VXn)*scrhR;)

                EXPAND(                                                       ,
                        usL_BXt = uL_BXt/rho*(uL_RHO*duL*duL - Bx*Bx)*scrhL; 
                        usR_BXt = uR_BXt/rho*(uR_RHO*duR*duR - Bx*Bx)*scrhR;  ,
                        usL_BXb = uL_BXb/rho*(uL_RHO*duL*duL - Bx*Bx)*scrhL;           
                        usR_BXb = uR_BXb/rho*(uR_RHO*duR*duR - Bx*Bx)*scrhR;)           

                if (S1L >= ZERO_F) {       //  ----  Region L*  ----

                    EXPAND(                                               ,
                    Flux(MXt,k,j,i) = fluxL_MXt + SL*(usL_MXt - uL_MXt);  ,
                    Flux(MXb,k,j,i) = fluxL_MXb + SL*(usL_MXb - uL_MXb);  
                    ) 
                    EXPAND(                                               ,
                    Flux(BXt,k,j,i) = fluxL_BXt + SL*(usL_BXt - uL_BXt);  ,
                    Flux(BXb,k,j,i) = fluxL_BXb + SL*(usL_BXb - uL_BXb);  
                    ) 

                }
                else if (S1R <= ZERO_F) {    //  ----  Region R*  ----
                
                    EXPAND(                                               ,
                    Flux(MXt,k,j,i) = fluxR_MXt + SR*(usR_MXt - uR_MXt);  ,
                    Flux(MXb,k,j,i) = fluxR_MXb + SR*(usR_MXb - uR_MXb);  
                    ) 
                    EXPAND(                                               ,
                    Flux(BXt,k,j,i) = fluxR_BXt + SR*(usR_BXt - uR_BXt);  ,
                    Flux(BXb,k,j,i) = fluxR_BXb + SR*(usR_BXb - uR_BXb);  
                    ) 
                    
                }
                else {

                    EXPAND(  ,
                        real usc_MXt;
                        real usc_BXt; ,
                        real usc_MXb;
                        real usc_BXb; )
                    
                //  Compute U** = Uc

                    sBx = (Bx > ZERO_F ? ONE_F : -ONE_F);

                    EXPAND(    ,
                        usc_MXt = HALF_F*(usR_MXt + usL_MXt + (usR_BXt - usL_BXt)*sBx*sqrho);  ,     
                        usc_MXb = HALF_F*(usR_MXb + usL_MXb + (usR_BXb - usL_BXb)*sBx*sqrho);  )
                    
                    EXPAND(    ,
                        usc_BXt = HALF_F*(   usR_BXt + usL_BXt + (usR_MXt - usL_MXt)*sBx/sqrho);  ,
                        usc_BXb = HALF_F*(   usR_BXb + usL_BXb + (usR_MXb - usL_MXb)*sBx/sqrho);  )

                    EXPAND(                                         ,
                        Flux(MXt,k,j,i) = usc_MXt*SM - Bx*usc_BXt;  ,
                        Flux(MXb,k,j,i) = usc_MXb*SM - Bx*usc_BXb;  )

                        
                    EXPAND(                                             ,
                        Flux(BXt,k,j,i) = usc_BXt*SM - Bx*usc_MXt/rho;  ,
                        Flux(BXb,k,j,i) = usc_BXb*SM - Bx*usc_MXb/rho;  )
                            
                }
                
            }

#endif
        }
                
        //6-- Compute maximum dt for this sweep
        const int ig = ioffset*i + joffset*j + koffset*k;

        invDt(k,j,i) = FMAX(cmax/dx(ig),invDt(k,j,i));

        // 7-- Store the flux in the emf components
        D_EXPAND(Et(k,j,i) = st*Flux(BXt,k,j,i); ,
                                                    ,
                    Eb(k,j,i) = sb*Flux(BXb,k,j,i); )
        
#if EMF_AVERAGE == UCT_CONTACT
        int s = 0;
        if (Flux(RHO,k,j,i) >  eps_UCT_CONTACT) s =  1;
        if (Flux(RHO,k,j,i) < -eps_UCT_CONTACT) s = -1;

        SV(k,j,i) = s;
#endif

    });


    idfx::popRegion();

}
