#include "../idefix.hpp"
#include "solversMHD.hpp"

#define ROE_AVERAGE 0

enum KWAVES {
    KFASTM = 0, KFASTP
#if HAVE_ENERGY
  , KENTRP
#endif

#if COMPONENTS >= 2
  , KSLOWM, KSLOWP
    #if COMPONENTS == 3
  , KALFVM, KALFVP
    #endif
#endif
};

/*! Return the sign of x. */
#define DSIGN(x)      ( (x) >= 0.0 ? (1.0) : (-1.0))

// Compute Riemann fluxes from states using ROE solver
void RoeMHD(DataBlock & data, int dir, real gamma, real C2Iso) {

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
    real delta    = 1.e-6;

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
        real dV[NVAR];

        // Conservative variables
        real uL[NVAR];
        real uR[NVAR];
        real dU[NVAR];

        // Flux (left and right)
        real fluxL[NVAR];
        real fluxR[NVAR];
        
        // Roe
        real Rc[NVAR][NVAR];
        real um[NVAR];

        // 1-- Store the primitive variables on the left, right, and averaged states
        for(int nv = 0 ; nv < NVAR; nv++) {
            vL[nv] = PrimL(nv,k,j,i);
            vR[nv] = PrimR(nv,k,j,i);
            dV[nv] = vR[nv] - vL[nv];
        }

        // 2-- Compute the conservative variables
        K_PrimToCons(uL, vL, gamma_m1);
        K_PrimToCons(uR, vR, gamma_m1);
        
        // 3-- Compute the left and right fluxes
        for(int nv = 0 ; nv < NVAR; nv++) {
            fluxL[nv] = uL[nv];
            fluxR[nv] = uR[nv];
            dU[nv] = uR[nv] - uL[nv];
        }
        
        K_Flux(fluxL, vL, fluxL, C2Iso, VXn, VXt, VXb, BXn, BXt, BXb, MXn);
        K_Flux(fluxR, vR, fluxR, C2Iso, VXn, VXt, VXb, BXn, BXt, BXb, MXn);
        
        // --- Compute the square of the sound speed
        real a, a2, a2L, a2R;
#if HAVE_ENERGY
        a2L = gamma * vL[PRS] / vL[RHO];
        a2R = gamma * vR[PRS] / vR[RHO];
        
#else
        a2L = C2Iso;
        a2R = C2Iso;
#endif
        
        // 5. Set eigenvectors components Rc = 0 initially  
        for(int nv1 = 0 ; nv1 < NVAR; nv1++) {
            for(int nv2 = 0 ; nv2 < NVAR; nv2++) {
                Rc[nv1][nv2] = 0;
            }
        }
        
        real sqr_rho_L, sqr_rho_R, sl, sr, rho, sqrt_rho;
        
        // 6c. Compute Roe averages 
        sqr_rho_L = sqrt(vL[RHO]);
        sqr_rho_R = sqrt(vR[RHO]);

        sl = sqr_rho_L/(sqr_rho_L + sqr_rho_R);
        sr = sqr_rho_R/(sqr_rho_L + sqr_rho_R);

        // sl = sr = 0.5;
        
        rho = sr*vL[RHO] + sl*vR[RHO];

        sqrt_rho = sqrt(rho);

        real u, v, w, Bx, By, Bz, sBx, bx, by, bz, bt2, b2, Btmag;
        
        EXPAND (u = sl*vL[VXn] + sr*vR[VXn];  ,
                v = sl*vL[VXt] + sr*vR[VXt];  ,
                w = sl*vL[VXb] + sr*vR[VXb];)

        EXPAND (Bx = sr*vL[BXn] + sl*vR[BXn];  ,
                By = sr*vL[BXt] + sl*vR[BXt];  ,
                Bz = sr*vL[BXb] + sl*vR[BXb];)

        EXPAND (Bx = sr*vL[BXn] + sl*vR[BXn];  ,
                By = sr*vL[BXt] + sl*vR[BXt];  ,
                Bz = sr*vL[BXb] + sl*vR[BXb];)

        sBx = (Bx >= 0.0 ? 1.0 : -1.0);

        EXPAND(bx = Bx/sqrt_rho;  ,
            by = By/sqrt_rho;  ,
            bz = Bz/sqrt_rho; )
        
        bt2   = EXPAND(0.0  , + by*by, + bz*bz);
        b2    = bx*bx + bt2;
        Btmag = sqrt(bt2*rho);

        real X, vdm, BdB;
        
        X  = EXPAND(dV[BXn]*dV[BXn], + dV[BXt]*dV[BXt], + dV[BXb]*dV[BXb]);
        X /= (sqr_rho_L + sqr_rho_R)*(sqr_rho_L + sqr_rho_R)*2.0;   

        vdm = EXPAND(u*dU[MXn],  + v*dU[MXt],  + w*dU[MXb]);

        BdB = EXPAND(Bx*dU[BXn], + By*dU[BXt], + Bz*dU[BXb]);

        
        real Bmag2L, Bmag2R, pL, pR;
        Bmag2L = EXPAND(vL[BX1]*vL[BX1] , + vL[BX2]*vL[BX2], + vL[BX3]*vL[BX3]);
        Bmag2R = EXPAND(vR[BX1]*vR[BX1] , + vR[BX2]*vR[BX2], + vR[BX3]*vR[BX3]);
#if HAVE_ENERGY
        pL  = vL[PRS] + HALF_F*Bmag2L;
        pR  = vR[PRS] + HALF_F*Bmag2R;
#else
        pL  = a2L*vL[RHO] + HALF_F*Bmag2L;
        pR  = a2R*vR[RHO] + HALF_F*Bmag2R;
#endif
        
        // 6d. Compute enthalpy and sound speed.
#if HAVE_ENERGY
        real vel2, HL, HR, H, Hgas;
        vel2    = EXPAND(u*u, + v*v, + w*w);
        dV[PRS] = gamma_m1*((0.5*vel2 - X)*dV[RHO] - vdm + dU[ENG] - BdB); 
        
        HL   = (uL[ENG] + pL)/vL[RHO];
        HR   = (uR[ENG] + pR)/vR[RHO];
        H    = sl*HL + sr*HR;   // total enthalpy

        Hgas = H - b2;         // gas enthalpy

        a2 = (2.0 - gamma)*X + gamma_m1*(Hgas - 0.5*vel2);
        if (a2 < 0.0) {
            //IDEFIX_ERROR("! Roe_Solver(): a2 < 0.0 !! \n");
        }
#else
        // in most cases a2L = a2R for isothermal MHD
        a2 = 0.5*(a2L + a2R) + X;
#endif
        
    /* ------------------------------------------------------------
        6e. Compute fast and slow magnetosonic speeds.

        The following expression appearing in the definitions
        of the fast magnetosonic speed 
        
        (a^2 - b^2)^2 + 4*a^2*bt^2 = (a^2 + b^2)^2 - 4*a^2*bx^2

        is always positive and avoids round-off errors.
        
        Note that we always use the total field to compute the 
        characteristic speeds.
        ------------------------------------------------------------ */
        
        real scrh, ca, cf, cs, ca2, cf2, cs2, alpha_f, alpha_s, beta_y, beta_z;
        scrh = a2 - b2;
        ca2  = bx*bx;
        scrh = scrh*scrh + 4.0*bt2*a2;    
        scrh = sqrt(scrh);    

        cf2 = 0.5*(a2 + b2 + scrh); 
        cs2 = a2*ca2/cf2;   // -- same as 0.5*(a2 + b2 - scrh)
        
        cf = sqrt(cf2);
        cs = sqrt(cs2);
        ca = sqrt(ca2);
        a  = sqrt(a2); 
        
        if (cf == cs) {
            alpha_f = 1.0;
            alpha_s = 0.0;
        }
        else if (a <= cs) {
            alpha_f = 0.0;
            alpha_s = 1.0;
        }
        else if (cf <= a) {
            alpha_f = 1.0;
            alpha_s = 0.0;
        }
        else{
            scrh    = 1.0/(cf2 - cs2);
            alpha_f = (a2  - cs2)*scrh;
            alpha_s = (cf2 -  a2)*scrh;
            alpha_f = FMAX(0.0, alpha_f);
            alpha_s = FMAX(0.0, alpha_s);
            alpha_f = sqrt(alpha_f);
            alpha_s = sqrt(alpha_s);
        }

        if (Btmag > 1.e-9) {
            SELECT(                  , 
                beta_y = DSIGN(By);  ,
                beta_y = By/Btmag; 
                beta_z = Bz/Btmag;)
        }
        else {
            SELECT(                        , 
                beta_y = 1.0;              ,
                beta_z = beta_y = 1.0;)
        }

    /* -------------------------------------------------------------------
        6f. Compute non-zero entries of conservative eigenvectors (Rc), 
            wave strength L*dU (=eta) for all 8 (or 7) waves using the
            expressions given by Eq. [4.18]--[4.21]. 
            Fast and slow eigenvectors are multiplied by a^2 while
            jumps are divided by a^2.
        
            Notes:
            - the expression on the paper has a typo in the very last term 
            of the energy component: it should be + and not - !
            - with background field splitting: additional terms must be 
            added to the energy component for fast, slow and Alfven waves.
            To obtain energy element, conservative eigenvector (with 
            total field) must be multiplied by | 0 0 0 0 -B0y -B0z 1 |.
            Also, H - b2 does not give gas enthalpy. A term b0*btot must 
            be added and eta (wave strength) should contain total field 
            and deviation's delta.
        ------------------------------------------------------------------- */

        // Fast wave:  u - c_f
        real lambda[NFLX], alambda[NFLX], eta[NFLX];
        real beta_dv, beta_dB, beta_v;
        
        int kk = KFASTM;
        lambda[kk] = u - cf;

        scrh    = alpha_s*cs*sBx;
        beta_dv = EXPAND(0.0, + beta_y*dV[VXt], + beta_z*dV[VXb]);
        beta_dB = EXPAND(0.0, + beta_y*dV[BXt], + beta_z*dV[BXb]);
        beta_v  = EXPAND(0.0, + beta_y*v,       + beta_z*w);

        Rc[RHO][kk] = alpha_f;
        EXPAND(Rc[MXn][kk] = alpha_f*lambda[kk];       ,
            Rc[MXt][kk] = alpha_f*v + scrh*beta_y;    ,
            Rc[MXb][kk] = alpha_f*w + scrh*beta_z;) 
        EXPAND(                                      ,
            Rc[BXt][kk] = alpha_s*a*beta_y/sqrt_rho;  ,
            Rc[BXb][kk] = alpha_s*a*beta_z/sqrt_rho;)

#if HAVE_ENERGY
        Rc[ENG][kk] =   alpha_f*(Hgas - u*cf) + scrh*beta_v
                    + alpha_s*a*Btmag/sqrt_rho;

        eta[kk] =   alpha_f*(X*dV[RHO] + dV[PRS]) + rho*scrh*beta_dv
                - rho*alpha_f*cf*dV[VXn]        + sqrt_rho*alpha_s*a*beta_dB;
#else
        eta[kk] =   alpha_f*(0.0*X + a2)*dV[RHO] + rho*scrh*beta_dv
                - rho*alpha_f*cf*dV[VXn] + sqrt_rho*alpha_s*a*beta_dB;
#endif
        
        eta[kk] *= 0.5/a2;

        // Fast wave:  u + c_f
        
        kk = KFASTP;
        lambda[kk] = u + cf;

        Rc[RHO][kk] = alpha_f;
        EXPAND( Rc[MXn][kk] = alpha_f*lambda[kk];        ,
                Rc[MXt][kk] = alpha_f*v - scrh*beta_y;  ,
                Rc[MXb][kk] = alpha_f*w - scrh*beta_z; ) 
        EXPAND(                                ,                                
                Rc[BXt][kk] = Rc[BXt][KFASTM];  ,
                Rc[BXb][kk] = Rc[BXb][KFASTM]; )

#if HAVE_ENERGY
        Rc[ENG][kk] =   alpha_f*(Hgas + u*cf) - scrh*beta_v
                    + alpha_s*a*Btmag/sqrt_rho;

        eta[kk] =   alpha_f*(X*dV[RHO] + dV[PRS]) - rho*scrh*beta_dv
                + rho*alpha_f*cf*dV[VXn]        + sqrt_rho*alpha_s*a*beta_dB;
#else
        eta[kk] =   alpha_f*(0.*X + a2)*dV[RHO] - rho*scrh*beta_dv
                + rho*alpha_f*cf*dV[VXn]      + sqrt_rho*alpha_s*a*beta_dB;
#endif

        eta[kk] *= 0.5/a2;

        // Entropy wave:  u
        
#if HAVE_ENERGY
        kk = KENTRP;
        lambda[kk] = u;

        Rc[RHO][kk] = 1.0;
        EXPAND( Rc[MXn][kk] = u; ,
                Rc[MXt][kk] = v; ,
                Rc[MXb][kk] = w; )
        Rc[ENG][kk] = 0.5*vel2 + (gamma - 2.0)/gamma_m1*X;

        eta[kk] = ((a2 - X)*dV[RHO] - dV[PRS])/a2;
#endif

    /* -----------------------------------------------------------------
        div.B wave (u): this wave exists when: 

        1) 8 wave formulation
        2) CT, since we always have 8 components, but it 
            carries zero jump.

        With GLM, KDIVB is replaced by KPSI_GLMM, KPSI_GLMP and these
        two waves should not enter in the Riemann solver (eta = 0.0) 
        since the 2x2 linear system formed by (B,psi) has already 
        been solved.
        ----------------------------------------------------------------- */

        /*kk = KDIVB;
        lambda[kk] = u;

        Rc[BXn][kk] = eta[kk] = 0.0;*/
        
#if COMPONENTS > 1    

        // Slow wave:  u - c_s
        
        scrh = alpha_f*cf*sBx;
        
        kk = KSLOWM;
        lambda[kk] = u - cs;

        Rc[RHO][kk] = alpha_s;
        EXPAND( Rc[MXn][kk] = alpha_s*lambda[kk];        ,
                Rc[MXt][kk] = alpha_s*v - scrh*beta_y;  ,
                Rc[MXb][kk] = alpha_s*w - scrh*beta_z;) 
        EXPAND(                                            ,                                
                Rc[BXt][kk] = - alpha_f*a*beta_y/sqrt_rho;  ,
                Rc[BXb][kk] = - alpha_f*a*beta_z/sqrt_rho; )

    #if HAVE_ENERGY
        Rc[ENG][kk] =   alpha_s*(Hgas - u*cs) - scrh*beta_v
                    - alpha_f*a*Btmag/sqrt_rho; 

        eta[kk] =   alpha_s*(X*dV[RHO] + dV[PRS]) - rho*scrh*beta_dv
                - rho*alpha_s*cs*dV[VXn]        - sqrt_rho*alpha_f*a*beta_dB;
    #else
        eta[kk] =   alpha_s*(0.*X + a2)*dV[RHO] - rho*scrh*beta_dv
                - rho*alpha_s*cs*dV[VXn]      - sqrt_rho*alpha_f*a*beta_dB;
    #endif

        eta[kk] *= 0.5/a2;

        // Slow wave:  u + c_s
        
        kk = KSLOWP;
        lambda[kk] = u + cs; 

        Rc[RHO][kk] = alpha_s;
        EXPAND(Rc[MXn][kk] = alpha_s*lambda[kk];         ,
                Rc[MXt][kk] = alpha_s*v + scrh*beta_y;   ,
                Rc[MXb][kk] = alpha_s*w + scrh*beta_z; ) 
        EXPAND(                                , 
                Rc[BXt][kk] = Rc[BXt][KSLOWM];   ,
                Rc[BXb][kk] = Rc[BXb][KSLOWM];)

    #if HAVE_ENERGY
        Rc[ENG][kk] =   alpha_s*(Hgas + u*cs) + scrh*beta_v
                    - alpha_f*a*Btmag/sqrt_rho;

        eta[kk] =   alpha_s*(X*dV[RHO] + dV[PRS]) + rho*scrh*beta_dv
                + rho*alpha_s*cs*dV[VXn]        - sqrt_rho*alpha_f*a*beta_dB; 
    #else
        eta[kk] =   alpha_s*(0.*X + a2)*dV[RHO] + rho*scrh*beta_dv
                + rho*alpha_s*cs*dV[VXn]      - sqrt_rho*alpha_f*a*beta_dB; 
    #endif

        eta[kk] *= 0.5/a2;

#endif // COMPONENTS > 1

#if COMPONENTS == 3

        // Alfven wave:  u - c_a

        kk = KALFVM;
        lambda[kk] = u - ca;

        Rc[MXt][kk] = - rho*beta_z;  
        Rc[MXb][kk] = + rho*beta_y;
        Rc[BXt][kk] = - sBx*sqrt_rho*beta_z;   
        Rc[BXb][kk] =   sBx*sqrt_rho*beta_y;
    #if HAVE_ENERGY
        Rc[ENG][kk] = - rho*(v*beta_z - w*beta_y);
    #endif

        eta[kk] = + beta_y*dV[VXb]               - beta_z*dV[VXt] 
                + sBx/sqrt_rho*(beta_y*dV[BXb] - beta_z*dV[BXt]);

        eta[kk] *= 0.5;

        // Alfven wave:  u + c_a 

        kk = KALFVP;
        lambda[kk] = u + ca;

        Rc[MXt][kk] = - Rc[MXt][KALFVM];  
        Rc[MXb][kk] = - Rc[MXb][KALFVM];
        Rc[BXt][kk] =   Rc[BXt][KALFVM];   
        Rc[BXb][kk] =   Rc[BXb][KALFVM];
    #if HAVE_ENERGY
        Rc[ENG][kk] = - Rc[ENG][KALFVM];
    #endif

        eta[kk] = - beta_y*dV[VXb]               + beta_z*dV[VXt] 
                + sBx/sqrt_rho*(beta_y*dV[BXb] - beta_z*dV[BXt]);

        eta[kk] *= 0.5;
#endif // COMPONENTS == 3

        // 6g. Compute maximum signal velocity

        real cmax = fabs (u) + cf;
        
        for(int nv = 0 ; nv < NVAR; nv++) {
            alambda[nv] = fabs(lambda[nv]);
        }

        // 6h. Entropy Fix 
        
        if (alambda[KFASTM] < 0.5*delta) {
            alambda[KFASTM] = lambda[KFASTM]*lambda[KFASTM]/delta + 0.25*delta;
        }
        if (alambda[KFASTP] < 0.5*delta) {
            alambda[KFASTP] = lambda[KFASTP]*lambda[KFASTP]/delta + 0.25*delta;
        }
#if COMPONENTS > 1
        if (alambda[KSLOWM] < 0.5*delta) {
            alambda[KSLOWM] = lambda[KSLOWM]*lambda[KSLOWM]/delta + 0.25*delta;
        }
        if (alambda[KSLOWP] < 0.5*delta) {
            alambda[KSLOWP] = lambda[KSLOWP]*lambda[KSLOWP]/delta + 0.25*delta; 
        }
#endif
    
        // 6i. Compute Roe numerical flux 
        
        for(int nv1 = 0 ; nv1 < NFLX; nv1++) {
            scrh = 0.0;
            for(int nv2 = 0 ; nv2 < NFLX; nv2++) {
                scrh += alambda[nv2]*eta[nv2]*Rc[nv1][nv2];
            }
            Flux(nv1,k,j,i) = 0.5*(fluxL[nv1] + fluxR[nv1] - scrh);
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
