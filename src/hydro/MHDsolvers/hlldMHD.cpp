#include "../idefix.hpp"
#include "solversMHD.hpp"

// Compute Riemann fluxes from states using HLLD solver
void HlldMHD(DataBlock & data, int dir, real gamma, real C2Iso) {

    int ioffset,joffset,koffset;
    int iextend, jextend,kextend;

    Kokkos::Profiling::pushRegion("TVDLF_Solver");
    
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


    real gamma_m1=gamma-ONE_F;

    // Define normal, tangent and bi-tanget indices

    int nDIR, tDIR, bDIR;
    int VXn, VXt, VXb;
    int BXn, BXt, BXb;
    int MXn, MXt, MXb;

    real st, sb;      // st and sb will be useful only when Hall is included
    switch(dir) {
        case(IDIR):
            ioffset = 1;
            D_EXPAND(               ,
                   jextend = 1;     ,
                   kextend = 1; )

            EXPAND(VXn = MXn = VX1; 
                   BXn = BX1;       , 
                   VXt = MXt = VX2; 
                   BXt = BX2;       , 
                   VXb = MXb = VX3;
                   BXb = BX3;       )

            Et = data.emf.ezi;
            Eb = data.emf.eyi;

            st = -1.0;
            sb = +1.0;
            break;
        case(JDIR):
            joffset=1;
            D_EXPAND( iextend = 1;  ,
                                    ,
                    kextend = 1;)
            EXPAND(VXn = MXn = VX2; 
                   BXn = BX2;       , 
                   VXt = MXt = VX1; 
                   BXt = BX1;       , 
                   VXb = MXb = VX3;
                   BXb = BX3;       )

            Et = data.emf.ezj;
            Eb = data.emf.exj;

            st = +1.0;
            sb = -1.0;
            break;
        case(KDIR):
            koffset=1;
            D_EXPAND( iextend = 1;  ,
                    jextend = 1;    ,
                    )
            EXPAND(VXn = MXn = VX3; 
                   BXn = BX3;       , 
                   VXt = MXt = VX1; 
                   BXt = BX1;       , 
                   VXb = MXb = VX2;
                   BXb = BX2;       )

            Et = data.emf.eyk;
            Eb = data.emf.exk;

            st = -1.0;
            sb = +1.0;
            break;
        default:
            IDEFIX_ERROR("Wrong direction");
    }

    nDIR = VXn-VX1; tDIR = VXt-VX1; bDIR = VXb-VX1;


    idefix_for("CalcRiemannFlux",data.beg[KDIR]-kextend,data.end[KDIR]+koffset+kextend,data.beg[JDIR]-jextend,data.end[JDIR]+joffset+jextend,data.beg[IDIR]-iextend,data.end[IDIR]+ioffset+iextend,
                KOKKOS_LAMBDA (int k, int j, int i) 
    {
        // Primitive variables
        real vL[NVAR];
        real vR[NVAR];

        // Conservative variables
        real uL[NVAR];
        real uR[NVAR];

        // Flux (left and right)
        real fluxL[NVAR];
        real fluxR[NVAR];

        // Signal speeds
        real cL, cR, cmax;

        // 1-- Store the primitive variables on the left, right, and averaged states
        for(int nv = 0 ; nv < NVAR; nv++) {
            vL[nv] = PrimL(nv,k,j,i);
            vR[nv] = PrimR(nv,k,j,i);
        }

        // 2-- Get the wave speed
        real gpr, b1, b2, b3, Btmag2, Bmag2;
#if HAVE_ENERGY
        gpr = gamma*vL[PRS];
#else
        gpr = C2Iso*vL[RHO];
#endif

        // -- get total field
        b1 = b2 = b3 = 0.0;
        EXPAND (b1 = vL[BXn];  ,
                b2 = vL[BXt];  ,
                b3 = vL[BXb];)

        Btmag2 = b2*b2 + b3*b3;
        Bmag2  = b1*b1 + Btmag2;

        cL = gpr - Bmag2;
        cL = gpr + Bmag2 + sqrt(cL*cL + 4.0*gpr*Btmag2);
        cL = sqrt(HALF_F*cL/vL[RHO]);
        
#if HAVE_ENERGY
        gpr = gamma*vR[PRS];
#else
        gpr = C2Iso*vR[RHO];
#endif

        // -- get total field
        b1 = b2 = b3 = 0.0;
        EXPAND (b1 = vR[BXn];  ,
                b2 = vR[BXt];  ,
                b3 = vR[BXb];)

        Btmag2 = b2*b2 + b3*b3;
        Bmag2  = b1*b1 + Btmag2;

        cR = gpr - Bmag2;
        cR = gpr + Bmag2 + sqrt(cR*cR + 4.0*gpr*Btmag2);
        cR = sqrt(HALF_F*cR/vR[RHO]);
        
        // 4.1 
        real cminL = vL[VXn] - cL;
        real cmaxL = vL[VXn] + cL;
        
        real cminR = vR[VXn] - cR;
        real cmaxR = vR[VXn] + cR;
        
        real SL = FMIN(cminL, cminR);
        real SR = FMAX(cmaxL, cmaxR);
        
        cmax  = FMAX(FABS(SL), FABS(SR));
        
        // 2-- Compute the conservative variables
        K_PrimToCons(uL, vL, gamma_m1);
        K_PrimToCons(uR, vR, gamma_m1);

        // 3-- Compute the left and right fluxes
        for(int nv = 0 ; nv < NVAR; nv++) {
            fluxL[nv] = uL[nv];
            fluxR[nv] = uR[nv];
        }
        
        K_Flux(fluxL, vL, fluxL, C2Iso, VXn, VXt, VXb, BXn, BXt, BXb, MXn);
        K_Flux(fluxR, vR, fluxR, C2Iso, VXn, VXt, VXb, BXn, BXt, BXb, MXn);
        
        real ptR, ptL;

#if HAVE_ENERGY
        ptL  = vL[PRS] + 0.5* ( EXPAND(vL[BX1]*vL[BX1] , + vL[BX2]*vL[BX2], + vL[BX3]*vL[BX3]) );
        ptR  = vR[PRS] + 0.5* ( EXPAND(vR[BX1]*vR[BX1] , + vR[BX2]*vR[BX2], + vR[BX3]*vR[BX3]) );
#else
        ptL  = C2Iso*vL[RHO] + 0.5* ( EXPAND(vL[BX1]*vL[BX1] , + vL[BX2]*vL[BX2], + vL[BX3]*vL[BX3]) );
        ptR  = C2Iso*vR[RHO] + 0.5* ( EXPAND(vR[BX1]*vR[BX1] , + vR[BX2]*vR[BX2], + vR[BX3]*vR[BX3]) );
#endif

        // 5-- Compute the flux from the left and right states
        if (SL > 0) {
            for (int nv = 0 ; nv < NFLX; nv++) {
                Flux(nv,k,j,i) = fluxL[nv];
            }
        }
        else if (SR < 0) {
            for (int nv = 0 ; nv < NFLX; nv++) {
                Flux(nv,k,j,i) = fluxR[nv];
            }
        }
        else {
            
            real usL[NVAR];
            real usR[NVAR];
            
            real scrh, scrhL, scrhR, duL, duR, sBx, Bx, Bx1, SM, S1L, S1R;
            
#if HAVE_ENERGY
            
            real Uhll[NVAR];
            real vs, pts, sqrL, sqrR, vsL, vsR, wsL, wsR;
            int revert_to_hllc;
            
            // 3c. Compute U*(L), U^*(R)

            scrh = 1.0/(SR - SL);
            Bx1  = Bx = (SR*vR[BXn] - SL*vL[BXn])*scrh; 
            sBx  = (Bx > 0.0 ? 1.0 : -1.0);

            duL  = SL - vL[VXn];
            duR  = SR - vR[VXn];

            scrh = 1.0/(duR*uR[RHO] - duL*uL[RHO]);
            SM   = (duR*uR[MXn] - duL*uL[MXn] - ptR + ptL)*scrh;

            pts  = duR*uR[RHO]*ptL - duL*uL[RHO]*ptR + 
                    vL[RHO]*vR[RHO]*duR*duL*(vR[VXn]- vL[VXn]);
            pts *= scrh;

            usL[RHO] = uL[RHO]*duL/(SL - SM);
            usR[RHO] = uR[RHO]*duR/(SR - SM);

            sqrL = sqrt(usL[RHO]);
            sqrR = sqrt(usR[RHO]);

            S1L = SM - fabs(Bx)/sqrL;
            S1R = SM + fabs(Bx)/sqrR;

            /* -----------------------------------------------------------------
            3d When S1L -> SL or S1R -> SR a degeneracy occurs. 
            Although Miyoshi & Kusano say that no jump exists, we don't
            think this is actually true. 
            Indeed, vy*, vz*, By*, Bz* cannot be solved independently. 
            In this case we revert to the HLLC solver of Li (2005),  except
            for the term v.B in the region, which we compute in our own way.
            Note, that by comparing the expressions of Li (2005) and
            Miyoshi & Kusano (2005), the only change involves a 
            re-definition of By* and Bz* in terms of By(HLL), Bz(HLL).
            ----------------------------------------------------------------- */

            revert_to_hllc = 0;

            if ( (S1L - SL) <  1.e-4*(SM - SL) ) revert_to_hllc = 1;
            if ( (S1R - SR) > -1.e-4*(SR - SM) ) revert_to_hllc = 1;
               
            if (revert_to_hllc) {

                scrh = 1.0/(SR - SL);
                
                for(int nv = 0 ; nv < NVAR; nv++) {
                    Uhll[nv]  = SR*uR[nv] - SL*uL[nv] + fluxL[nv] - fluxR[nv];
                    Uhll[nv] *= scrh;
                }
                
                // WHERE'S THE PRESSURE ?!?!?!?
                EXPAND(usL[BXn] = usR[BXn] = Uhll[BXn];   ,
                    usL[BXt] = usR[BXt] = Uhll[BXt];   ,
                    usL[BXb] = usR[BXb] = Uhll[BXb];)
        
                S1L = S1R = SM; // region ** should never be computed since
                                // fluxes are given in terms of UL* and UR*

            }
            else {

            // 3e. Compute states in the * regions

                scrhL = (uL[RHO]*duL*duL - Bx*Bx)/(uL[RHO]*duL*(SL - SM) - Bx*Bx);
                scrhR = (uR[RHO]*duR*duR - Bx*Bx)/(uR[RHO]*duR*(SR - SM) - Bx*Bx);
        
                EXPAND(usL[BXn]  = Bx1;         ,
                    usL[BXt]  = uL[BXt]*scrhL;  ,
                    usL[BXb]  = uL[BXb]*scrhL;)           

                EXPAND(usR[BXn] = Bx1;          ,
                    usR[BXt] = uR[BXt]*scrhR;   ,   
                    usR[BXb] = uR[BXb]*scrhR;)     
                
            }

            scrhL = Bx/(uL[RHO]*duL);
            scrhR = Bx/(uR[RHO]*duR);

            EXPAND(                                          ;  ,
                    vsL = vL[VXt] - scrhL*(usL[BXt] - uL[BXt]);
                    vsR = vR[VXt] - scrhR*(usR[BXt] - uR[BXt]);  ,

                    wsL = vL[VXb] - scrhL*(usL[BXb] - uL[BXb]);
                    wsR = vR[VXb] - scrhR*(usR[BXb] - uR[BXb]); )

            EXPAND(usL[MXn] = usL[RHO]*SM; 
                    usR[MXn] = usR[RHO]*SM;   ,
            
                    usL[MXt] = usL[RHO]*vsL;
                    usR[MXt] = usR[RHO]*vsR;  ,

                    usL[MXb] = usL[RHO]*wsL;
                    usR[MXb] = usR[RHO]*wsR;)

            /* -- Energy -- */

            scrhL  = EXPAND(vL[VXn]*Bx1, + vL[VXt]*uL[BXt], + vL[VXb]*uL[BXb]);
            scrhL -= EXPAND(     SM*Bx1, +    vsL*usL[BXt], +    wsL*usL[BXb]);
            usL[ENG]  = duL*uL[ENG] - ptL*vL[VXn] + pts*SM + Bx*scrhL;
            usL[ENG] /= SL - SM;

            scrhR  = EXPAND(vR[VXn]*Bx1, + vR[VXt]*uR[BXt], + vR[VXb]*uR[BXb]);
            scrhR -= EXPAND(     SM*Bx1, +    vsR*usR[BXt], +    wsR*usR[BXb]);
            usR[ENG] = duR*uR[ENG] - ptR*vR[VXn] + pts*SM + Bx*scrhR;
            usR[ENG] /= SR - SM;

        // 3c. Compute flux when S1L > 0 or S1R < 0

            if (S1L >= 0.0) {       //  ----  Region L*

                for(int nv = 0 ; nv < NVAR; nv++) {
                    Flux(nv,k,j,i) = fluxL[nv] + SL*(usL[nv] - uL[nv]);
                }

            }else if (S1R <= 0.0) {    //  ----  Region R*
            
                for(int nv = 0 ; nv < NVAR; nv++) {
                    Flux(nv,k,j,i) = fluxR[nv] + SR*(usR[nv] - uR[nv]);
                }
                
            } else {   // -- This state exists only if B_x != 0

                // Compute U**
                real vss, wss;
                real ussl[NVAR];
                real ussr[NVAR];

                ussl[RHO] = usL[RHO];
                ussr[RHO] = usR[RHO];
                
                EXPAND(                           ,
            
                    vss  = sqrL*vsL + sqrR*vsR + (usR[BXt] - usL[BXt])*sBx;       
                    vss /= sqrL + sqrR;        ,
                    
                    wss  = sqrL*wsL + sqrR*wsR + (usR[BXb] - usL[BXb])*sBx;
                    wss /= sqrL + sqrR;)

                EXPAND(ussl[MXn] = ussl[RHO]*SM;
                    ussr[MXn] = ussr[RHO]*SM;    ,
            
                    ussl[MXt] = ussl[RHO]*vss;
                    ussr[MXt] = ussr[RHO]*vss;  ,
                
                    ussl[MXb] = ussl[RHO]*wss;
                    ussr[MXb] = ussr[RHO]*wss;)           
            
                EXPAND(ussl[BXn] = ussr[BXn] = Bx1;   ,

                    ussl[BXt]  = sqrL*usR[BXt] + sqrR*usL[BXt] + sqrL*sqrR*(vsR - vsL)*sBx;
                    ussl[BXt] /= sqrL + sqrR;        
                    ussr[BXt]  = ussl[BXt];        ,
                
                    ussl[BXb]  = sqrL*usR[BXb] + sqrR*usL[BXb] + sqrL*sqrR*(wsR - wsL)*sBx;
                    ussl[BXb] /= sqrL + sqrR;        
                    ussr[BXb]  = ussl[BXb];)
                
                // -- Energy jump

                scrhL  = EXPAND(SM*Bx1, +  vsL*usL [BXt], +  wsL*usL [BXb]);
                scrhL -= EXPAND(SM*Bx1, +  vss*ussl[BXt], +  wss*ussl[BXb]);

                scrhR  = EXPAND(SM*Bx1, +  vsR*usR [BXt], +  wsR*usR [BXb]);
                scrhR -= EXPAND(SM*Bx1, +  vss*ussr[BXt], +  wss*ussr[BXb]);

                ussl[ENG] = usL[ENG] - sqrL*scrhL*sBx;
                ussr[ENG] = usR[ENG] + sqrR*scrhR*sBx;


                if (SM >= 0.0) {           //  ----  Region L**
                    for(int nv = 0 ; nv < NVAR; nv++) {
                        Flux(nv,k,j,i) = fluxL[nv] + S1L*(ussl[nv]  - usL[nv])
                                                    + SL*(usL[nv] - uL[nv]);
                    }
                }else{                   //  ----  Region R**
                    for(int nv = 0 ; nv < NVAR; nv++) {
                        Flux(nv,k,j,i) = fluxR[nv] + S1R*(ussr[nv]  - usR[nv])
                                                    + SR*(usR[nv] - uR[nv]);
                    }
                }
            }  // end if (S1L < 0 S1R > 0)
#else
            real usc[NVAR];
            real rho, sqrho;
            int revert_to_hll;
            
            scrh = 1.0/(SR - SL);
            duL = SL - vL[VXn];
            duR = SR - vR[VXn];

            Bx1 = Bx = (SR*vR[BXn] - SL*vL[BXn])*scrh; 

            rho                = (uR[RHO]*duR - uL[RHO]*duL)*scrh;
            Flux(RHO,k,j,i) = (SL*uR[RHO]*duR - SR*uL[RHO]*duL)*scrh;
                
        /* ---------------------------
                compute S*
            --------------------------- */

            sqrho = sqrt(rho);

            SM  = Flux(RHO,k,j,i)/rho;
            S1L = SM - fabs(Bx)/sqrho;
            S1R = SM + fabs(Bx)/sqrho;

            /* ---------------------------------------------
                Prevent degeneracies when S1L -> SL or 
                S1R -> SR. Revert to HLL if necessary.
            --------------------------------------------- */

            revert_to_hll = 0;

            if ( (S1L - SL) <  1.e-4*(SR - SL) ) revert_to_hll = 1;
            if ( (S1R - SR) > -1.e-4*(SR - SL) ) revert_to_hll = 1;

            if (revert_to_hll) {
                scrh = 1.0/(SR - SL);
                for(int nv = 0 ; nv < NVAR; nv++) {
                    Flux(nv,k,j,i) = SL*SR*(uR[nv] - uL[nv]) +
                                        SR*fluxL[nv] - SL*fluxR[nv];
                    Flux(nv,k,j,i) *= scrh;
                }
            }
            else {

                Flux(MXn,k,j,i) = (SR*fluxL[MXn] - SL*fluxR[MXn] 
                                        + SR*SL*(uR[MXn] - uL[MXn]))*scrh;

                Flux(BXn,k,j,i) = SR*SL*(uR[BXn] - uL[BXn])*scrh;


            /* ---------------------------
                        Compute U*  
                --------------------------- */
                
                scrhL = 1.0/((SL - S1L)*(SL - S1R));
                scrhR = 1.0/((SR - S1L)*(SR - S1R));

                EXPAND(                                                        ;  ,
                        usL[MXt] = rho*vL[VXt] - Bx*uL[BXt]*(SM - vL[VXn])*scrhL;
                        usR[MXt] = rho*vR[VXt] - Bx*uR[BXt]*(SM - vR[VXn])*scrhR;  ,

                        usL[MXb] = rho*vL[VXb] - Bx*uL[BXb]*(SM - vL[VXn])*scrhL;
                        usR[MXb] = rho*vR[VXb] - Bx*uR[BXb]*(SM - vR[VXn])*scrhR;)

                EXPAND(                                                      ;  ,
                        usL[BXt] = uL[BXt]/rho*(uL[RHO]*duL*duL - Bx*Bx)*scrhL; 
                        usR[BXt] = uR[BXt]/rho*(uR[RHO]*duR*duR - Bx*Bx)*scrhR;  ,

                        usL[BXb] = uL[BXb]/rho*(uL[RHO]*duL*duL - Bx*Bx)*scrhL;           
                        usR[BXb] = uR[BXb]/rho*(uR[RHO]*duR*duR - Bx*Bx)*scrhR;)           

                if (S1L >= 0.0) {       /*  ----  Region L*  ---- */

                    EXPAND(                                                    ;  ,
                    Flux(MXt,k,j,i) = fluxL[MXt] + SL*(usL[MXt] - uL[MXt]);  ,
                    Flux(MXb,k,j,i) = fluxL[MXb] + SL*(usL[MXb] - uL[MXb]);  
                    ) 
                    EXPAND(                                                    ;  ,
                    Flux(BXt,k,j,i) = fluxL[BXt] + SL*(usL[BXt] - uL[BXt]);  ,
                    Flux(BXb,k,j,i) = fluxL[BXb] + SL*(usL[BXb] - uL[BXb]);  
                    ) 

                }
                else if (S1R <= 0.0) {    /*  ----  Region R*  ---- */
                
                    EXPAND(                                                    ;  ,
                    Flux(MXt,k,j,i) = fluxR[MXt] + SR*(usR[MXt] - uR[MXt]);  ,
                    Flux(MXb,k,j,i) = fluxR[MXb] + SR*(usR[MXb] - uR[MXb]);  
                    ) 
                    EXPAND(                                                    ;  ,
                    Flux(BXt,k,j,i) = fluxR[BXt] + SR*(usR[BXt] - uR[BXt]);  ,
                    Flux(BXb,k,j,i) = fluxR[BXb] + SR*(usR[BXb] - uR[BXb]);  
                    ) 
                    
                } else {
                                
                /* ---------------------------
                        Compute U** = Uc
                    --------------------------- */

                    sBx = (Bx > 0.0 ? 1.0 : -1.0);

                    EXPAND(                                                  ;  ,
                        usc[MXt] = 0.5*(usR[MXt] + usL[MXt] 
                                        + (usR[BXt] - usL[BXt])*sBx*sqrho);  ,     
                        usc[MXb] = 0.5*(   usR[MXb] + usL[MXb] 
                                        + (usR[BXb] - usL[BXb])*sBx*sqrho);)
                    
                    EXPAND(                                                  ;  ,
                        usc[BXt] = 0.5*(   usR[BXt] + usL[BXt]  
                                        + (usR[MXt] - usL[MXt])*sBx/sqrho);  ,
                        usc[BXb] = 0.5*(   usR[BXb] + usL[BXb] 
                                        + (usR[MXb] - usL[MXb])*sBx/sqrho);)

                    EXPAND(                                               ;  ,
                        Flux(MXt,k,j,i) = usc[MXt]*SM - Bx*usc[BXt];  ,
                        Flux(MXb,k,j,i) = usc[MXb]*SM - Bx*usc[BXb]; )

                        
                    EXPAND(                                                   ;  ,
                        Flux(BXt,k,j,i) = usc[BXt]*SM - Bx*usc[MXt]/rho;  ,
                        Flux(BXb,k,j,i) = usc[BXb]*SM - Bx*usc[MXb]/rho;)
                            
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

    });


    Kokkos::Profiling::popRegion();

}
