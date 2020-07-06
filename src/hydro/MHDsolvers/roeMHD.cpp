#include "../idefix.hpp"
#include "solversMHD.hpp"

#define ROE_AVERAGE 0


/*! Return the sign of x. */
#define DSIGN(x)      ( (x) >= 0.0 ? (1.0) : (-1.0))

// Compute Riemann fluxes from states using ROE solver
void RoeMHD(DataBlock & data, int dir, real gamma, real C2Iso) {

    int ioffset,joffset,koffset;
    int iextend, jextend,kextend;

    Kokkos::Profiling::pushRegion("ROE_MHD");
    
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
    real delta    = 1.e-6;

    // Define normal, tangent and bi-tanget indices

    EXPAND( int BXt;  ,
                      ,
            int BXb;  )

    // st and sb will be useful only when Hall is included
    D_EXPAND( real st;  ,
                        ,
              real sb;  )
    
    switch(dir) {
        case(IDIR):
            ioffset = 1;
            D_EXPAND(               ,
                   jextend = 1;     ,
                   kextend = 1; )

            EXPAND(/*BXn = BX1;*/   , 
                   BXt = BX2;       , 
                   BXb = BX3;       )

            Et = data.emf.ezi;
            Eb = data.emf.eyi;
            SV = data.emf.svx;

            D_EXPAND( st = -ONE_F;  ,
                                    ,
                      sb = +ONE_F;  )
            break;
        case(JDIR):
            joffset=1;
            D_EXPAND( iextend = 1;  ,
                                    ,
                    kextend = 1;)
            EXPAND(/*BXn = BX2;*/   , 
                   BXt = BX1;       , 
                   BXb = BX3;       )

            Et = data.emf.ezj;
            Eb = data.emf.exj;
            SV = data.emf.svy;

            D_EXPAND( st = +ONE_F;  ,
                                    ,
                      sb = -ONE_F;  )
            break;
        case(KDIR):
            koffset=1;
            D_EXPAND( iextend = 1;  ,
                    jextend = 1;    ,
                    )
            EXPAND(/*BXn = BX3;*/   , 
                   BXt = BX1;       , 
                   BXb = BX2;       )

            Et = data.emf.eyk;
            Eb = data.emf.exk;
            SV = data.emf.svz;

            D_EXPAND( st = -ONE_F;  ,
                                    ,
                      sb = +ONE_F;  )
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
        
        // --- Compute the square of the sound speed
        real a, a2;
#if HAVE_ENERGY
        // real a2L = gamma * vL_PRS / vL_RHO;
        // real a2R = gamma * vR_PRS / vR_RHO;
        
#else
        real a2L = C2Iso;
        real a2R = C2Iso;
#endif

        real dV_RHO = vR_RHO - vL_RHO;
        EXPAND(
            real dV_VX1 = vR_VX1 - vL_VX1;
            real dV_BX1 = vR_BX1 - vL_BX1; ,
            real dV_VX2 = vR_VX2 - vL_VX2;
            real dV_BX2 = vR_BX2 - vL_BX2; ,
            real dV_VX3 = vR_VX3 - vL_VX3;
            real dV_BX3 = vR_BX3 - vL_BX3; )

#if HAVE_ENERGY
        // real dU_RHO = uR_RHO - uL_RHO;
        EXPAND(
            real dU_MX1 = uR_MX1 - uL_MX1;
            real dU_BX1 = uR_BX1 - uL_BX1; ,
            real dU_MX2 = uR_MX2 - uL_MX2;
            real dU_BX2 = uR_BX2 - uL_BX2; ,
            real dU_MX3 = uR_MX3 - uL_MX3;
            real dU_BX3 = uR_BX3 - uL_BX3; )
        
        real dV_PRS = vR_PRS - vL_PRS;
        real dU_ENG = uR_ENG - uL_ENG;
#endif

        // 5. Set eigenvectors components Rc = 0 initially  
        real Rc_RHO_KFASTM = ZERO_F;
        EXPAND(
            real Rc_RHO_KFASTP = ZERO_F;
            real Rc_RHO_KDIVB = ZERO_F; ,
            real Rc_RHO_KSLOWM = ZERO_F;
            real Rc_RHO_KSLOWP = ZERO_F; ,
            real Rc_RHO_KALFVM = ZERO_F;
            real Rc_RHO_KALFVP = ZERO_F; )
        
        EXPAND(
            real Rc_MX1_KFASTM = ZERO_F;
            EXPAND(
                real Rc_MX1_KFASTP = ZERO_F;
                real Rc_MX1_KDIVB = ZERO_F; ,
                real Rc_MX1_KSLOWM = ZERO_F;
                real Rc_MX1_KSLOWP = ZERO_F; ,
                real Rc_MX1_KALFVM = ZERO_F;
                real Rc_MX1_KALFVP = ZERO_F; )
            real Rc_BX1_KFASTM = ZERO_F;
            EXPAND(
                real Rc_BX1_KFASTP = ZERO_F;
                real Rc_BX1_KDIVB = ZERO_F; ,
                real Rc_BX1_KSLOWM = ZERO_F;
                real Rc_BX1_KSLOWP = ZERO_F; ,
                real Rc_BX1_KALFVM = ZERO_F;
                real Rc_BX1_KALFVP = ZERO_F; )
            ,
            real Rc_MX2_KFASTM = ZERO_F;
            EXPAND(
                real Rc_MX2_KFASTP = ZERO_F;
                real Rc_MX2_KDIVB = ZERO_F; ,
                real Rc_MX2_KSLOWM = ZERO_F;
                real Rc_MX2_KSLOWP = ZERO_F; ,
                real Rc_MX2_KALFVM = ZERO_F;
                real Rc_MX2_KALFVP = ZERO_F; )
            real Rc_BX2_KFASTM = ZERO_F;
            EXPAND(
                real Rc_BX2_KFASTP = ZERO_F;
                real Rc_BX2_KDIVB = ZERO_F; ,
                real Rc_BX2_KSLOWM = ZERO_F;
                real Rc_BX2_KSLOWP = ZERO_F; ,
                real Rc_BX2_KALFVM = ZERO_F;
                real Rc_BX2_KALFVP = ZERO_F; )
            ,
            real Rc_MX3_KFASTM = ZERO_F;
            EXPAND(
                real Rc_MX3_KFASTP = ZERO_F;
                real Rc_MX3_KDIVB = ZERO_F; ,
                real Rc_MX3_KSLOWM = ZERO_F;
                real Rc_MX3_KSLOWP = ZERO_F; ,
                real Rc_MX3_KALFVM = ZERO_F;
                real Rc_MX3_KALFVP = ZERO_F; )
            real Rc_BX3_KFASTM = ZERO_F;
            EXPAND(
                real Rc_BX3_KFASTP = ZERO_F;
                real Rc_BX3_KDIVB = ZERO_F; ,
                real Rc_BX3_KSLOWM = ZERO_F;
                real Rc_BX3_KSLOWP = ZERO_F; ,
                real Rc_BX3_KALFVM = ZERO_F;
                real Rc_BX3_KALFVP = ZERO_F; )
        )

#if HAVE_ENERGY
        real Rc_RHO_KENTRP = ZERO_F;
        EXPAND(
            real Rc_MX1_KENTRP = ZERO_F;
            real Rc_BX1_KENTRP = ZERO_F; ,
            real Rc_MX2_KENTRP = ZERO_F;
            real Rc_BX2_KENTRP = ZERO_F; ,
            real Rc_MX3_KENTRP = ZERO_F;
            real Rc_BX3_KENTRP = ZERO_F; )
        
        real Rc_ENG_KFASTM = ZERO_F;
        EXPAND(
            real Rc_ENG_KFASTP = ZERO_F;
            real Rc_ENG_KDIVB = ZERO_F; ,
            real Rc_ENG_KSLOWM = ZERO_F;
            real Rc_ENG_KSLOWP = ZERO_F; ,
            real Rc_ENG_KALFVM = ZERO_F;
            real Rc_ENG_KALFVP = ZERO_F; )
        real Rc_ENG_KENTRP = ZERO_F;
        
        EXPAND(
            real dU_MXn; real dU_BXn;  ,
            real dU_MXt; real dU_BXt;  ,
            real dU_MXb; real dU_BXb;  )
        
        EXPAND (
        if (dir == IDIR) {
            EXPAND (
            dU_MXn = dU_MX1;
            dU_BXn = dU_BX1;  ,
            dU_MXt = dU_MX2;
            dU_BXt = dU_BX2;  ,
            dU_MXb = dU_MX3;
            dU_BXb = dU_BX3;  )
        }
        ,
        if (dir == JDIR) {
            EXPAND (
            dU_MXn = dU_MX2;
            dU_BXn = dU_BX2;  ,
            dU_MXt = dU_MX1;
            dU_BXt = dU_BX1;  ,
            dU_MXb = dU_MX3;
            dU_BXb = dU_BX3;  )
        }
        ,
        if (dir == KDIR) {
            EXPAND (
            dU_MXn = dU_MX3;
            dU_BXn = dU_BX3;  ,
            dU_MXt = dU_MX1;
            dU_BXt = dU_BX1;  ,
            dU_MXb = dU_MX2;
            dU_BXb = dU_BX2;  )
        }
        )
        
#endif
        
        EXPAND(
            real vL_VXn; real vR_VXn;
            real vL_BXn; real vR_BXn;
            real dV_VXn; real dV_BXn;  ,
            real vL_VXt; real vR_VXt;
            real vL_BXt; real vR_BXt;
            real dV_VXt; real dV_BXt;  ,
            real vL_VXb; real vR_VXb;
            real vL_BXb; real vR_BXb;
            real dV_VXb; real dV_BXb;  )
        
        EXPAND (
        if (dir == IDIR) {
            EXPAND (
            vL_VXn = vL_VX1;
            vR_VXn = vR_VX1;
            vL_BXn = vL_BX1;
            vR_BXn = vR_BX1;
            dV_VXn = dV_VX1;
            dV_BXn = dV_BX1;  ,
            vL_VXt = vL_VX2;
            vR_VXt = vR_VX2;
            vL_BXt = vL_BX2;
            vR_BXt = vR_BX2;
            dV_VXt = dV_VX2;
            dV_BXt = dV_BX2;  ,
            vL_VXb = vL_VX3;
            vR_VXb = vR_VX3;
            vL_BXb = vL_BX3;
            vR_BXb = vR_BX3;
            dV_VXb = dV_VX3;
            dV_BXb = dV_BX3;  )
        }
        ,
        if (dir == JDIR) {
            EXPAND (
            vL_VXn = vL_VX2;
            vR_VXn = vR_VX2;
            vL_BXn = vL_BX2;
            vR_BXn = vR_BX2;
            dV_VXn = dV_VX2;
            dV_BXn = dV_BX2;  ,
            vL_VXt = vL_VX1;
            vR_VXt = vR_VX1;
            vL_BXt = vL_BX1;
            vR_BXt = vR_BX1;
            dV_VXt = dV_VX1;
            dV_BXt = dV_BX1;  ,
            vL_VXb = vL_VX3;
            vR_VXb = vR_VX3;
            vL_BXb = vL_BX3;
            vR_BXb = vR_BX3;
            dV_VXb = dV_VX3;
            dV_BXb = dV_BX3;  )
        }
        ,
        if (dir == KDIR) {
            EXPAND (
            vL_VXn = vL_VX3;
            vR_VXn = vR_VX3;
            vL_BXn = vL_BX3;
            vR_BXn = vR_BX3;
            dV_VXn = dV_VX3;
            dV_BXn = dV_BX3;  ,
            vL_VXt = vL_VX1;
            vR_VXt = vR_VX1;
            vL_BXt = vL_BX1;
            vR_BXt = vR_BX1;
            dV_VXt = dV_VX1;
            dV_BXt = dV_BX1;  ,
            vL_VXb = vL_VX2;
            vR_VXb = vR_VX2;
            vL_BXb = vL_BX2;
            vR_BXb = vR_BX2;
            dV_VXb = dV_VX2;
            dV_BXb = dV_BX2;  )
        }
        )
        
        real sqr_rho_L, sqr_rho_R, sl, sr, rho, sqrt_rho;
        
        // 6c. Compute Roe averages 
        sqr_rho_L = sqrt(vL_RHO);
        sqr_rho_R = sqrt(vR_RHO);

        sl = sqr_rho_L/(sqr_rho_L + sqr_rho_R);
        sr = sqr_rho_R/(sqr_rho_L + sqr_rho_R);

        // sl = sr = HALF_F;
        
        rho = sr*vL_RHO + sl*vR_RHO;

        sqrt_rho = sqrt(rho);

        real sBx, bt2, b2, Btmag;
        
        EXPAND (real u = sl*vL_VXn + sr*vR_VXn;  ,
                real v = sl*vL_VXt + sr*vR_VXt;  ,
                real w = sl*vL_VXb + sr*vR_VXb;  )

        EXPAND (real Bx = sr*vL_BXn + sl*vR_BXn;  ,
                real By = sr*vL_BXt + sl*vR_BXt;  ,
                real Bz = sr*vL_BXb + sl*vR_BXb;  )

        sBx = (Bx >= ZERO_F ? ONE_F : -ONE_F);

        EXPAND(real bx = Bx/sqrt_rho;  ,
            real by = By/sqrt_rho;     ,
            real bz = Bz/sqrt_rho;     )
        
        bt2   = EXPAND(ZERO_F  , + by*by, + bz*bz);
        b2    = bx*bx + bt2;
        Btmag = sqrt(bt2*rho);

        
        real X  = EXPAND(dV_BXn*dV_BXn, + dV_BXt*dV_BXt, + dV_BXb*dV_BXb);
        X /= (sqr_rho_L + sqr_rho_R)*(sqr_rho_L + sqr_rho_R)*2.0;   

#if HAVE_ENERGY
        real Bmag2L, Bmag2R;
        Bmag2L = EXPAND(vL_BX1*vL_BX1 , + vL_BX2*vL_BX2, + vL_BX3*vL_BX3);
        Bmag2R = EXPAND(vR_BX1*vR_BX1 , + vR_BX2*vR_BX2, + vR_BX3*vR_BX3);
//         
        real pL  = vL_PRS + HALF_F*Bmag2L;
        real pR  = vR_PRS + HALF_F*Bmag2R;
#else
        // real pL  = a2L*vL_RHO + HALF_F*Bmag2L;
        // real pR  = a2R*vR_RHO + HALF_F*Bmag2R;
#endif
        
        // 6d. Compute enthalpy and sound speed.
#if HAVE_ENERGY
        real vdm = EXPAND(u*dU_MXn,  + v*dU_MXt,  + w*dU_MXb);

        real BdB = EXPAND(Bx*dU_BXn, + By*dU_BXt, + Bz*dU_BXb);
        
        real vel2, HL, HR, H, Hgas;
        vel2    = EXPAND(u*u, + v*v, + w*w);
        dV_PRS = gamma_m1*((HALF_F*vel2 - X)*dV_RHO - vdm + dU_ENG - BdB); 
        
        HL   = (uL_ENG + pL)/vL_RHO;
        HR   = (uR_ENG + pR)/vR_RHO;
        H    = sl*HL + sr*HR;   // total enthalpy

        Hgas = H - b2;         // gas enthalpy

        a2 = (2.0 - gamma)*X + gamma_m1*(Hgas - HALF_F*vel2);
        if (a2 < ZERO_F) {
            //IDEFIX_ERROR("! Roe_Solver(): a2 < 0.0 !! \n");
        }
#else
        // in most cases a2L = a2R for isothermal MHD
        a2 = HALF_F*(a2L + a2R) + X;
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

        cf2 = HALF_F*(a2 + b2 + scrh); 
        cs2 = a2*ca2/cf2;   // -- same as 0.5*(a2 + b2 - scrh)
        
        cf = sqrt(cf2);
        cs = sqrt(cs2);
        ca = sqrt(ca2);
        a  = sqrt(a2); 
        
        if (cf == cs) {
            alpha_f = ONE_F;
            alpha_s = ZERO_F;
        }
        else if (a <= cs) {
            alpha_f = ZERO_F;
            alpha_s = ONE_F;
        }
        else if (cf <= a) {
            alpha_f = ONE_F;
            alpha_s = ZERO_F;
        }
        else{
            scrh    = ONE_F/(cf2 - cs2);
            alpha_f = (a2  - cs2)*scrh;
            alpha_s = (cf2 -  a2)*scrh;
            alpha_f = FMAX(ZERO_F, alpha_f);
            alpha_s = FMAX(ZERO_F, alpha_s);
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
            SELECT(                     , 
                beta_y = ONE_F;         ,
                beta_z = beta_y = ONE_F;)
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
        
        EXPAND(
            real Rc_MXn_KFASTM = ZERO_F;
            EXPAND(
                real Rc_MXn_KFASTP = ZERO_F;
                real Rc_MXn_KDIVB = ZERO_F;  ,
                real Rc_MXn_KSLOWM = ZERO_F;
                real Rc_MXn_KSLOWP = ZERO_F; ,
                real Rc_MXn_KALFVM = ZERO_F;
                real Rc_MXn_KALFVP = ZERO_F; )
            real Rc_BXn_KFASTM = ZERO_F;
            EXPAND(
                real Rc_BXn_KFASTP = ZERO_F;
                real Rc_BXn_KDIVB = ZERO_F;  ,
                real Rc_BXn_KSLOWM = ZERO_F;
                real Rc_BXn_KSLOWP = ZERO_F; ,
                real Rc_BXn_KALFVM = ZERO_F;
                real Rc_BXn_KALFVP = ZERO_F; )
            ,
            real Rc_MXt_KFASTM = ZERO_F;
            EXPAND(
                real Rc_MXt_KFASTP = ZERO_F;
                real Rc_MXt_KDIVB = ZERO_F;  ,
                real Rc_MXt_KSLOWM = ZERO_F;
                real Rc_MXt_KSLOWP = ZERO_F; ,
                real Rc_MXt_KALFVM = ZERO_F;
                real Rc_MXt_KALFVP = ZERO_F; )
            real Rc_BXt_KFASTM = ZERO_F;
            EXPAND(
                real Rc_BXt_KFASTP = ZERO_F;
                real Rc_BXt_KDIVB = ZERO_F;  ,
                real Rc_BXt_KSLOWM = ZERO_F;
                real Rc_BXt_KSLOWP = ZERO_F; ,
                real Rc_BXt_KALFVM = ZERO_F;
                real Rc_BXt_KALFVP = ZERO_F; )
            ,
            real Rc_MXb_KFASTM = ZERO_F;
            EXPAND(
                real Rc_MXb_KFASTP = ZERO_F;
                real Rc_MXb_KDIVB = ZERO_F;  ,
                real Rc_MXb_KSLOWM = ZERO_F;
                real Rc_MXb_KSLOWP = ZERO_F; ,
                real Rc_MXb_KALFVM = ZERO_F;
                real Rc_MXb_KALFVP = ZERO_F; )
            real Rc_BXb_KFASTM = ZERO_F;
            EXPAND(
                real Rc_BXb_KFASTP = ZERO_F;
                real Rc_BXb_KDIVB = ZERO_F;  ,
                real Rc_BXb_KSLOWM = ZERO_F;
                real Rc_BXb_KSLOWP = ZERO_F; ,
                real Rc_BXb_KALFVM = ZERO_F;
                real Rc_BXb_KALFVP = ZERO_F; )
        )

#if HAVE_ENERGY
        EXPAND(
            real Rc_MXn_KENTRP = ZERO_F;
            real Rc_BXn_KENTRP = ZERO_F; ,
            real Rc_MXt_KENTRP = ZERO_F;
            real Rc_BXt_KENTRP = ZERO_F; ,
            real Rc_MXb_KENTRP = ZERO_F;
            real Rc_BXb_KENTRP = ZERO_F; )
#endif

        real lambda_KFASTM;
        EXPAND(
            real lambda_KFASTP;
            real lambda_KDIVB;  ,
            real lambda_KSLOWM;
            real lambda_KSLOWP; ,
            real lambda_KALFVM;
            real lambda_KALFVP; )
        
        real eta_KFASTM;
        EXPAND(
            real eta_KFASTP;
            real eta_KDIVB;  ,
            real eta_KSLOWM;
            real eta_KSLOWP; ,
            real eta_KALFVM;
            real eta_KALFVP; )
        
#if HAVE_ENERGY
        real lambda_KENTRP;
        real eta_KENTRP;
#endif

        real beta_dv, beta_dB;
        
        // Fast wave:  u - c_f
        
        lambda_KFASTM = u - cf;

        scrh    = alpha_s*cs*sBx;
        beta_dv = EXPAND(ZERO_F, + beta_y*dV_VXt, + beta_z*dV_VXb);
        beta_dB = EXPAND(ZERO_F, + beta_y*dV_BXt, + beta_z*dV_BXb);

        Rc_RHO_KFASTM = alpha_f;
        EXPAND( Rc_MXn_KFASTM = alpha_f*lambda_KFASTM;    ,
                Rc_MXt_KFASTM = alpha_f*v + scrh*beta_y;  ,
                Rc_MXb_KFASTM = alpha_f*w + scrh*beta_z;  ) 
        EXPAND( Rc_BXt_KFASTM = ZERO_F;                     ,
                Rc_BXt_KFASTM = alpha_s*a*beta_y/sqrt_rho;  ,
                Rc_BXb_KFASTM = alpha_s*a*beta_z/sqrt_rho;  )

#if HAVE_ENERGY
        real beta_v  = EXPAND(ZERO_F, + beta_y*v,      + beta_z*w);
        
        Rc_ENG_KFASTM =   alpha_f*(Hgas - u*cf) + scrh*beta_v
                    + alpha_s*a*Btmag/sqrt_rho;

        eta_KFASTM = alpha_f*(X*dV_RHO + dV_PRS);
#else
        eta_KFASTM = alpha_f*a2*dV_RHO;
                
#endif
        eta_KFASTM += rho*scrh*beta_dv - rho*alpha_f*cf*dV_VXn + sqrt_rho*alpha_s*a*beta_dB;
        eta_KFASTM *= HALF_F/a2;

        // Fast wave:  u + c_f

        lambda_KFASTP = u + cf;

        Rc_RHO_KFASTP = alpha_f;
        EXPAND( Rc_MXn_KFASTP = alpha_f*lambda_KFASTP;    ,
                Rc_MXt_KFASTP = alpha_f*v - scrh*beta_y;  ,
                Rc_MXb_KFASTP = alpha_f*w - scrh*beta_z;  ) 
        EXPAND( Rc_BXn_KFASTP = ZERO_F;         ,                                
                Rc_BXt_KFASTP = Rc_BXt_KFASTM;  ,
                Rc_BXb_KFASTP = Rc_BXb_KFASTM;  )

#if HAVE_ENERGY
        Rc_ENG_KFASTP =   alpha_f*(Hgas + u*cf) - scrh*beta_v
                    + alpha_s*a*Btmag/sqrt_rho;

        eta_KFASTP =   alpha_f*(X*dV_RHO + dV_PRS);
#else
        eta_KFASTP =   alpha_f*a2*dV_RHO;
#endif
        eta_KFASTP += rho*alpha_f*cf*dV_VXn + sqrt_rho*alpha_s*a*beta_dB - rho*scrh*beta_dv;
        eta_KFASTP *= HALF_F/a2;

        // Entropy wave:  u
        
#if HAVE_ENERGY
        lambda_KENTRP = u;

        Rc_RHO_KENTRP = ONE_F;
        EXPAND( Rc_MXn_KENTRP = u; ,
                Rc_MXt_KENTRP = v; ,
                Rc_MXb_KENTRP = w; )
        EXPAND( Rc_BXn_KENTRP = ZERO_F; ,
                Rc_BXt_KENTRP = ZERO_F; ,
                Rc_BXb_KENTRP = ZERO_F; )
        Rc_ENG_KENTRP = HALF_F*vel2 + (gamma - 2.0)/gamma_m1*X;

        eta_KENTRP = ((a2 - X)*dV_RHO - dV_PRS)/a2;
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

        lambda_KDIVB = u;
        eta_KDIVB = ZERO_F;
        
#if COMPONENTS > 1    

        // Slow wave:  u - c_s
        
        scrh = alpha_f*cf*sBx;

        lambda_KSLOWM = u - cs;

        Rc_RHO_KSLOWM = alpha_s;
        EXPAND( Rc_MXn_KSLOWM = alpha_s*lambda_KSLOWM;    ,
                Rc_MXt_KSLOWM = alpha_s*v - scrh*beta_y;  ,
                Rc_MXb_KSLOWM = alpha_s*w - scrh*beta_z;  ) 
        EXPAND( Rc_BXn_KSLOWM = ZERO_F;                       ,                                
                Rc_BXt_KSLOWM = - alpha_f*a*beta_y/sqrt_rho;  ,
                Rc_BXb_KSLOWM = - alpha_f*a*beta_z/sqrt_rho;  )

    #if HAVE_ENERGY
        Rc_ENG_KSLOWM =   alpha_s*(Hgas - u*cs) - scrh*beta_v
                    - alpha_f*a*Btmag/sqrt_rho; 

        eta_KSLOWM = alpha_s*(X*dV_RHO + dV_PRS);
    #else
        eta_KSLOWM = alpha_s*a2*dV_RHO;
    #endif
        eta_KSLOWM += - rho*scrh*beta_dv - rho*alpha_s*cs*dV_VXn - sqrt_rho*alpha_f*a*beta_dB;
        eta_KSLOWM *= HALF_F/a2;

        // Slow wave:  u + c_s
        
        lambda_KSLOWP = u + cs; 

        Rc_RHO_KSLOWP = alpha_s;
        EXPAND( Rc_MXn_KSLOWP = alpha_s*lambda_KSLOWP;    ,
                Rc_MXt_KSLOWP = alpha_s*v + scrh*beta_y;  ,
                Rc_MXb_KSLOWP = alpha_s*w + scrh*beta_z;  ) 
        EXPAND( Rc_BXt_KSLOWP = ZERO_F;         , 
                Rc_BXt_KSLOWP = Rc_BXt_KSLOWM;  ,
                Rc_BXb_KSLOWP = Rc_BXb_KSLOWM;  )

    #if HAVE_ENERGY
        Rc_ENG_KSLOWP =   alpha_s*(Hgas + u*cs) + scrh*beta_v
                    - alpha_f*a*Btmag/sqrt_rho;

        eta_KSLOWP =   alpha_s*(X*dV_RHO + dV_PRS);
    #else
        eta_KSLOWP =   alpha_s*a2*dV_RHO;
                 
    #endif
        eta_KSLOWP += rho*scrh*beta_dv + rho*alpha_s*cs*dV_VXn - sqrt_rho*alpha_f*a*beta_dB;
        eta_KSLOWP *= HALF_F/a2;

#endif // COMPONENTS > 1

#if COMPONENTS == 3

        // Alfven wave:  u - c_a

        lambda_KALFVM = u - ca;

        Rc_MXt_KALFVM = - rho*beta_z;  
        Rc_MXb_KALFVM = + rho*beta_y;
        Rc_BXt_KALFVM = - sBx*sqrt_rho*beta_z;   
        Rc_BXb_KALFVM =   sBx*sqrt_rho*beta_y;
    #if HAVE_ENERGY
        Rc_ENG_KALFVM = - rho*(v*beta_z - w*beta_y);
    #endif

        eta_KALFVM = + beta_y*dV_VXb - beta_z*dV_VXt 
                + sBx/sqrt_rho*(beta_y*dV_BXb - beta_z*dV_BXt);

        eta_KALFVM *= HALF_F;

        // Alfven wave:  u + c_a 

        lambda_KALFVP = u + ca;

        Rc_MXt_KALFVP = - Rc_MXt_KALFVM;  
        Rc_MXb_KALFVP = - Rc_MXb_KALFVM;
        Rc_BXt_KALFVP =   Rc_BXt_KALFVM;   
        Rc_BXb_KALFVP =   Rc_BXb_KALFVM;
    #if HAVE_ENERGY
        Rc_ENG_KALFVP = - Rc_ENG_KALFVM;
    #endif

        eta_KALFVP = - beta_y*dV_VXb + beta_z*dV_VXt 
                + sBx/sqrt_rho*(beta_y*dV_BXb - beta_z*dV_BXt);

        eta_KALFVP *= HALF_F;
#endif // COMPONENTS == 3

        // 6g. Compute maximum signal velocity

        real cmax = FABS(u) + cf;
        
        real alambda_KFASTM = FABS(lambda_KFASTM);
        EXPAND(
            real alambda_KFASTP = FABS(lambda_KFASTP);
            real alambda_KDIVB  = FABS(lambda_KDIVB);   ,
            real alambda_KSLOWM = FABS(lambda_KSLOWM);
            real alambda_KSLOWP = FABS(lambda_KSLOWP); ,
            real alambda_KALFVM = FABS(lambda_KALFVM);
            real alambda_KALFVP = FABS(lambda_KALFVP); )
#if HAVE_ENERGY
        real alambda_KENTRP = FABS(lambda_KENTRP);
#endif

        // 6h. Entropy Fix 
        
        if (alambda_KFASTM < HALF_F*delta) {
            alambda_KFASTM = lambda_KFASTM*lambda_KFASTM/delta + 0.25*delta;
        }
        if (alambda_KFASTP < HALF_F*delta) {
            alambda_KFASTP = lambda_KFASTP*lambda_KFASTP/delta + 0.25*delta;
        }
#if COMPONENTS > 1
        if (alambda_KSLOWM < HALF_F*delta) {
            alambda_KSLOWM = lambda_KSLOWM*lambda_KSLOWM/delta + 0.25*delta;
        }
        if (alambda_KSLOWP < HALF_F*delta) {
            alambda_KSLOWP = lambda_KSLOWP*lambda_KSLOWP/delta + 0.25*delta; 
        }
#endif
    
        // 6i. Compute Roe numerical flux 
        
        // MXi
        EXPAND (
        if (dir == IDIR) {
            EXPAND (
            //MX1 = MXn;
                Rc_MX1_KFASTM = Rc_MXn_KFASTM;
                EXPAND(
                    Rc_MX1_KFASTP = Rc_MXn_KFASTP;
                    Rc_MX1_KDIVB  = Rc_MXn_KDIVB;  ,
                    Rc_MX1_KSLOWM = Rc_MXn_KSLOWM;
                    Rc_MX1_KSLOWP = Rc_MXn_KSLOWP; ,
                    Rc_MX1_KALFVM = Rc_MXn_KALFVM;
                    Rc_MX1_KALFVP = Rc_MXn_KALFVP; )
#if HAVE_ENERGY
                Rc_MX1_KENTRP = Rc_MXn_KENTRP;
#endif
              ,
            
            //MX2 = MXt;
                Rc_MX2_KFASTM = Rc_MXt_KFASTM;
                EXPAND(
                    Rc_MX2_KFASTP = Rc_MXt_KFASTP;
                    Rc_MX2_KDIVB  = Rc_MXt_KDIVB;  ,
                    Rc_MX2_KSLOWM = Rc_MXt_KSLOWM;
                    Rc_MX2_KSLOWP = Rc_MXt_KSLOWP; ,
                    Rc_MX2_KALFVM = Rc_MXt_KALFVM;
                    Rc_MX2_KALFVP = Rc_MXt_KALFVP; )
#if HAVE_ENERGY
                Rc_MX2_KENTRP = Rc_MXt_KENTRP;
#endif
              ,
            //MX3 = MXb;
                Rc_MX3_KFASTM = Rc_MXb_KFASTM;
                EXPAND(
                    Rc_MX3_KFASTP = Rc_MXb_KFASTP;
                    Rc_MX3_KDIVB  = Rc_MXb_KDIVB;  ,
                    Rc_MX3_KSLOWM = Rc_MXb_KSLOWM;
                    Rc_MX3_KSLOWP = Rc_MXb_KSLOWP; ,
                    Rc_MX3_KALFVM = Rc_MXb_KALFVM;
                    Rc_MX3_KALFVP = Rc_MXb_KALFVP; )
#if HAVE_ENERGY
                Rc_MX3_KENTRP = Rc_MXb_KENTRP;
#endif
            )
        }
        ,
        if (dir == JDIR) {
            EXPAND (
            //MX2 = MXn;
                Rc_MX2_KFASTM = Rc_MXn_KFASTM;
                EXPAND(
                    Rc_MX2_KFASTP = Rc_MXn_KFASTP;
                    Rc_MX2_KDIVB  = Rc_MXn_KDIVB;  ,
                    Rc_MX2_KSLOWM = Rc_MXn_KSLOWM;
                    Rc_MX2_KSLOWP = Rc_MXn_KSLOWP; ,
                    Rc_MX2_KALFVM = Rc_MXn_KALFVM;
                    Rc_MX2_KALFVP = Rc_MXn_KALFVP; )
#if HAVE_ENERGY
                Rc_MX2_KENTRP = Rc_MXn_KENTRP;
#endif
              ,
            
            //MX1 = MXt;
                Rc_MX1_KFASTM = Rc_MXt_KFASTM;
                EXPAND(
                    Rc_MX1_KFASTP = Rc_MXt_KFASTP;
                    Rc_MX1_KDIVB  = Rc_MXt_KDIVB;  ,
                    Rc_MX1_KSLOWM = Rc_MXt_KSLOWM;
                    Rc_MX1_KSLOWP = Rc_MXt_KSLOWP; ,
                    Rc_MX1_KALFVM = Rc_MXt_KALFVM;
                    Rc_MX1_KALFVP = Rc_MXt_KALFVP; )
#if HAVE_ENERGY
                Rc_MX1_KENTRP = Rc_MXt_KENTRP;
#endif
              ,
            //MX3 = MXb;
                Rc_MX3_KFASTM = Rc_MXb_KFASTM;
                EXPAND(
                    Rc_MX3_KFASTP = Rc_MXb_KFASTP;
                    Rc_MX3_KDIVB  = Rc_MXb_KDIVB;  ,
                    Rc_MX3_KSLOWM = Rc_MXb_KSLOWM;
                    Rc_MX3_KSLOWP = Rc_MXb_KSLOWP; ,
                    Rc_MX3_KALFVM = Rc_MXb_KALFVM;
                    Rc_MX3_KALFVP = Rc_MXb_KALFVP; )
#if HAVE_ENERGY
                Rc_MX3_KENTRP = Rc_MXb_KENTRP;
#endif
            )
        }
        ,
        if (dir == KDIR) {
            EXPAND (
            //MX3 = MXn;
                Rc_MX3_KFASTM = Rc_MXn_KFASTM;
                EXPAND(
                    Rc_MX3_KFASTP = Rc_MXn_KFASTP;
                    Rc_MX3_KDIVB  = Rc_MXn_KDIVB;  ,
                    Rc_MX3_KSLOWM = Rc_MXn_KSLOWM;
                    Rc_MX3_KSLOWP = Rc_MXn_KSLOWP; ,
                    Rc_MX3_KALFVM = Rc_MXn_KALFVM;
                    Rc_MX3_KALFVP = Rc_MXn_KALFVP; )
#if HAVE_ENERGY
                Rc_MX3_KENTRP = Rc_MXn_KENTRP;
#endif
              ,
            
            //MX1 = MXt;
                Rc_MX1_KFASTM = Rc_MXt_KFASTM;
                EXPAND(
                    Rc_MX1_KFASTP = Rc_MXt_KFASTP;
                    Rc_MX1_KDIVB  = Rc_MXt_KDIVB;  ,
                    Rc_MX1_KSLOWM = Rc_MXt_KSLOWM;
                    Rc_MX1_KSLOWP = Rc_MXt_KSLOWP; ,
                    Rc_MX1_KALFVM = Rc_MXt_KALFVM;
                    Rc_MX1_KALFVP = Rc_MXt_KALFVP; )
#if HAVE_ENERGY
                Rc_MX1_KENTRP = Rc_MXt_KENTRP;
#endif
              ,
            //MX2 = MXb;
                Rc_MX2_KFASTM = Rc_MXb_KFASTM;
                EXPAND(
                    Rc_MX2_KFASTP = Rc_MXb_KFASTP;
                    Rc_MX2_KDIVB  = Rc_MXb_KDIVB;  ,
                    Rc_MX2_KSLOWM = Rc_MXb_KSLOWM;
                    Rc_MX2_KSLOWP = Rc_MXb_KSLOWP; ,
                    Rc_MX2_KALFVM = Rc_MXb_KALFVM;
                    Rc_MX2_KALFVP = Rc_MXb_KALFVP; )
#if HAVE_ENERGY
                Rc_MX2_KENTRP = Rc_MXb_KENTRP;
#endif
            )
        }
        )
        
        // BXi
        EXPAND (
        if (dir == IDIR) {
            EXPAND (
            //BX1 = BXn;
                Rc_BX1_KFASTM = Rc_BXn_KFASTM;
                EXPAND(
                    Rc_BX1_KFASTP = Rc_BXn_KFASTP;
                    Rc_BX1_KDIVB  = Rc_BXn_KDIVB;  ,
                    Rc_BX1_KSLOWM = Rc_BXn_KSLOWM;
                    Rc_BX1_KSLOWP = Rc_BXn_KSLOWP; ,
                    Rc_BX1_KALFVM = Rc_BXn_KALFVM;
                    Rc_BX1_KALFVP = Rc_BXn_KALFVP; )
#if HAVE_ENERGY
                Rc_BX1_KENTRP = Rc_BXn_KENTRP;
#endif
              ,
            
            //BX2 = BXt;
                Rc_BX2_KFASTM = Rc_BXt_KFASTM;
                EXPAND(
                    Rc_BX2_KFASTP = Rc_BXt_KFASTP;
                    Rc_BX2_KDIVB  = Rc_BXt_KDIVB;  ,
                    Rc_BX2_KSLOWM = Rc_BXt_KSLOWM;
                    Rc_BX2_KSLOWP = Rc_BXt_KSLOWP; ,
                    Rc_BX2_KALFVM = Rc_BXt_KALFVM;
                    Rc_BX2_KALFVP = Rc_BXt_KALFVP; )
#if HAVE_ENERGY
                Rc_BX2_KENTRP = Rc_BXt_KENTRP;
#endif
              ,
            //BX3 = BXb;
                Rc_BX3_KFASTM = Rc_BXb_KFASTM;
                EXPAND(
                    Rc_BX3_KFASTP = Rc_BXb_KFASTP;
                    Rc_BX3_KDIVB  = Rc_BXb_KDIVB;  ,
                    Rc_BX3_KSLOWM = Rc_BXb_KSLOWM;
                    Rc_BX3_KSLOWP = Rc_BXb_KSLOWP; ,
                    Rc_BX3_KALFVM = Rc_BXb_KALFVM;
                    Rc_BX3_KALFVP = Rc_BXb_KALFVP; )
#if HAVE_ENERGY
                Rc_BX3_KENTRP = Rc_BXb_KENTRP;
#endif
            )
        }
        ,
        if (dir == JDIR) {
            EXPAND (
            //BX2 = BXn;
                Rc_BX2_KFASTM = Rc_BXn_KFASTM;
                EXPAND(
                    Rc_BX2_KFASTP = Rc_BXn_KFASTP;
                    Rc_BX2_KDIVB  = Rc_BXn_KDIVB;  ,
                    Rc_BX2_KSLOWM = Rc_BXn_KSLOWM;
                    Rc_BX2_KSLOWP = Rc_BXn_KSLOWP; ,
                    Rc_BX2_KALFVM = Rc_BXn_KALFVM;
                    Rc_BX2_KALFVP = Rc_BXn_KALFVP; )
#if HAVE_ENERGY
                Rc_BX2_KENTRP = Rc_BXn_KENTRP;
#endif
              ,
            
            //BX1 = BXt;
                Rc_BX1_KFASTM = Rc_BXt_KFASTM;
                EXPAND(
                    Rc_BX1_KFASTP = Rc_BXt_KFASTP;
                    Rc_BX1_KDIVB  = Rc_BXt_KDIVB;  ,
                    Rc_BX1_KSLOWM = Rc_BXt_KSLOWM;
                    Rc_BX1_KSLOWP = Rc_BXt_KSLOWP; ,
                    Rc_BX1_KALFVM = Rc_BXt_KALFVM;
                    Rc_BX1_KALFVP = Rc_BXt_KALFVP; )
#if HAVE_ENERGY
                Rc_BX1_KENTRP = Rc_BXt_KENTRP;
#endif
              ,
            //BX3 = BXb;
                Rc_BX3_KFASTM = Rc_BXb_KFASTM;
                EXPAND(
                    Rc_BX3_KFASTP = Rc_BXb_KFASTP;
                    Rc_BX3_KDIVB  = Rc_BXb_KDIVB;  ,
                    Rc_BX3_KSLOWM = Rc_BXb_KSLOWM;
                    Rc_BX3_KSLOWP = Rc_BXb_KSLOWP; ,
                    Rc_BX3_KALFVM = Rc_BXb_KALFVM;
                    Rc_BX3_KALFVP = Rc_BXb_KALFVP; )
#if HAVE_ENERGY
                Rc_BX3_KENTRP = Rc_BXb_KENTRP;
#endif
            )
        }
        ,
        if (dir == KDIR) {
            EXPAND (
            //BX3 = BXn;
                Rc_BX3_KFASTM = Rc_BXn_KFASTM;
                EXPAND(
                    Rc_BX3_KFASTP = Rc_BXn_KFASTP;
                    Rc_BX3_KDIVB  = Rc_BXn_KDIVB;  ,
                    Rc_BX3_KSLOWM = Rc_BXn_KSLOWM;
                    Rc_BX3_KSLOWP = Rc_BXn_KSLOWP; ,
                    Rc_BX3_KALFVM = Rc_BXn_KALFVM;
                    Rc_BX3_KALFVP = Rc_BXn_KALFVP; )
#if HAVE_ENERGY
                Rc_BX3_KENTRP = Rc_BXn_KENTRP;
#endif
              ,
            
            //BX1 = BXt;
                Rc_BX1_KFASTM = Rc_BXt_KFASTM;
                EXPAND(
                    Rc_BX1_KFASTP = Rc_BXt_KFASTP;
                    Rc_BX1_KDIVB  = Rc_BXt_KDIVB;  ,
                    Rc_BX1_KSLOWM = Rc_BXt_KSLOWM;
                    Rc_BX1_KSLOWP = Rc_BXt_KSLOWP; ,
                    Rc_BX1_KALFVM = Rc_BXt_KALFVM;
                    Rc_BX1_KALFVP = Rc_BXt_KALFVP; )
#if HAVE_ENERGY
                Rc_BX1_KENTRP = Rc_BXt_KENTRP;
#endif
              ,
            //BX2 = BXb;
                Rc_BX2_KFASTM = Rc_BXb_KFASTM;
                EXPAND(
                    Rc_BX2_KFASTP = Rc_BXb_KFASTP;
                    Rc_BX2_KDIVB  = Rc_BXb_KDIVB;  ,
                    Rc_BX2_KSLOWM = Rc_BXb_KSLOWM;
                    Rc_BX2_KSLOWP = Rc_BXb_KSLOWP; ,
                    Rc_BX2_KALFVM = Rc_BXb_KALFVM;
                    Rc_BX2_KALFVP = Rc_BXb_KALFVP; )
#if HAVE_ENERGY
                Rc_BX2_KENTRP = Rc_BXb_KENTRP;
#endif
            )
        }
        )
            
            scrh = alambda_KFASTM*eta_KFASTM*Rc_RHO_KFASTM;
            EXPAND(
                scrh += alambda_KFASTP*eta_KFASTP*Rc_RHO_KFASTP;
                scrh += alambda_KDIVB *eta_KDIVB *Rc_RHO_KDIVB;   ,
                scrh += alambda_KSLOWM*eta_KSLOWM*Rc_RHO_KSLOWM; 
                scrh += alambda_KSLOWP*eta_KSLOWP*Rc_RHO_KSLOWP;  ,
                scrh += alambda_KALFVM*eta_KALFVM*Rc_RHO_KALFVM; 
                scrh += alambda_KALFVP*eta_KALFVP*Rc_RHO_KALFVP;  )
#if HAVE_ENERGY
            scrh += alambda_KENTRP*eta_KENTRP*Rc_RHO_KENTRP; 
#endif
            
            Flux(RHO,k,j,i) = HALF_F*(fluxL_RHO + fluxR_RHO - scrh);
            
        EXPAND ( scrh = alambda_KFASTM*eta_KFASTM*Rc_MX1_KFASTM;
                 EXPAND(
                    scrh += alambda_KFASTP*eta_KFASTP*Rc_MX1_KFASTP;
                    scrh += alambda_KDIVB *eta_KDIVB *Rc_MX1_KDIVB;   ,
                    scrh += alambda_KSLOWM*eta_KSLOWM*Rc_MX1_KSLOWM; 
                    scrh += alambda_KSLOWP*eta_KSLOWP*Rc_MX1_KSLOWP;  ,
                    scrh += alambda_KALFVM*eta_KALFVM*Rc_MX1_KALFVM; 
                    scrh += alambda_KALFVP*eta_KALFVP*Rc_MX1_KALFVP;  )
#if HAVE_ENERGY
                 scrh += alambda_KENTRP*eta_KENTRP*Rc_MX1_KENTRP; 
#endif
                 Flux(MX1,k,j,i) = HALF_F*(fluxL_MX1 + fluxR_MX1 - scrh);
                 
                 scrh = alambda_KFASTM*eta_KFASTM*Rc_BX1_KFASTM;
                 EXPAND(
                    scrh += alambda_KFASTP*eta_KFASTP*Rc_BX1_KFASTP;
                    scrh += alambda_KDIVB *eta_KDIVB *Rc_BX1_KDIVB;   ,
                    scrh += alambda_KSLOWM*eta_KSLOWM*Rc_BX1_KSLOWM; 
                    scrh += alambda_KSLOWP*eta_KSLOWP*Rc_BX1_KSLOWP;  ,
                    scrh += alambda_KALFVM*eta_KALFVM*Rc_BX1_KALFVM; 
                    scrh += alambda_KALFVP*eta_KALFVP*Rc_BX1_KALFVP;  )
#if HAVE_ENERGY
                 scrh += alambda_KENTRP*eta_KENTRP*Rc_BX1_KENTRP; 
#endif
                 Flux(BX1,k,j,i) = HALF_F*(fluxL_BX1 + fluxR_BX1 - scrh);    ,
                 
                 scrh = alambda_KFASTM*eta_KFASTM*Rc_MX2_KFASTM;
                 EXPAND(
                    scrh += alambda_KFASTP*eta_KFASTP*Rc_MX2_KFASTP;
                    scrh += alambda_KDIVB *eta_KDIVB *Rc_MX2_KDIVB;   ,
                    scrh += alambda_KSLOWM*eta_KSLOWM*Rc_MX2_KSLOWM; 
                    scrh += alambda_KSLOWP*eta_KSLOWP*Rc_MX2_KSLOWP;  ,
                    scrh += alambda_KALFVM*eta_KALFVM*Rc_MX2_KALFVM; 
                    scrh += alambda_KALFVP*eta_KALFVP*Rc_MX2_KALFVP;  )
#if HAVE_ENERGY
                 scrh += alambda_KENTRP*eta_KENTRP*Rc_MX2_KENTRP; 
#endif
                 Flux(MX2,k,j,i) = HALF_F*(fluxL_MX2 + fluxR_MX2 - scrh);
                 
                 scrh = alambda_KFASTM*eta_KFASTM*Rc_BX2_KFASTM;
                 EXPAND(
                    scrh += alambda_KFASTP*eta_KFASTP*Rc_BX2_KFASTP;
                    scrh += alambda_KDIVB *eta_KDIVB *Rc_BX2_KDIVB;   ,
                    scrh += alambda_KSLOWM*eta_KSLOWM*Rc_BX2_KSLOWM; 
                    scrh += alambda_KSLOWP*eta_KSLOWP*Rc_BX2_KSLOWP;  ,
                    scrh += alambda_KALFVM*eta_KALFVM*Rc_BX2_KALFVM; 
                    scrh += alambda_KALFVP*eta_KALFVP*Rc_BX2_KALFVP;  )
#if HAVE_ENERGY
                 scrh += alambda_KENTRP*eta_KENTRP*Rc_BX2_KENTRP; 
#endif
                 Flux(BX2,k,j,i) = HALF_F*(fluxL_BX2 + fluxR_BX2 - scrh);    ,
                 
                 scrh = alambda_KFASTM*eta_KFASTM*Rc_MX3_KFASTM;
                 EXPAND(
                    scrh += alambda_KFASTP*eta_KFASTP*Rc_MX3_KFASTP;
                    scrh += alambda_KDIVB *eta_KDIVB *Rc_MX3_KDIVB;   ,
                    scrh += alambda_KSLOWM*eta_KSLOWM*Rc_MX3_KSLOWM; 
                    scrh += alambda_KSLOWP*eta_KSLOWP*Rc_MX3_KSLOWP;  ,
                    scrh += alambda_KALFVM*eta_KALFVM*Rc_MX3_KALFVM; 
                    scrh += alambda_KALFVP*eta_KALFVP*Rc_MX3_KALFVP;  )
#if HAVE_ENERGY
                 scrh += alambda_KENTRP*eta_KENTRP*Rc_MX3_KENTRP; 
#endif
                 Flux(MX3,k,j,i) = HALF_F*(fluxL_MX3 + fluxR_MX3 - scrh);
                 
                 scrh = alambda_KFASTM*eta_KFASTM*Rc_BX3_KFASTM;
                 EXPAND(
                    scrh += alambda_KFASTP*eta_KFASTP*Rc_BX3_KFASTP;
                    scrh += alambda_KDIVB *eta_KDIVB *Rc_BX3_KDIVB;   ,
                    scrh += alambda_KSLOWM*eta_KSLOWM*Rc_BX3_KSLOWM; 
                    scrh += alambda_KSLOWP*eta_KSLOWP*Rc_BX3_KSLOWP;  ,
                    scrh += alambda_KALFVM*eta_KALFVM*Rc_BX3_KALFVM; 
                    scrh += alambda_KALFVP*eta_KALFVP*Rc_BX3_KALFVP;  )
#if HAVE_ENERGY
                 scrh += alambda_KENTRP*eta_KENTRP*Rc_BX3_KENTRP; 
#endif
                 Flux(BX3,k,j,i) = HALF_F*(fluxL_BX3 + fluxR_BX3 - scrh);    )

#if HAVE_ENERGY
        scrh = alambda_KFASTM*eta_KFASTM*Rc_ENG_KFASTM;
        EXPAND(
            scrh += alambda_KFASTP*eta_KFASTP*Rc_ENG_KFASTP;
            scrh += alambda_KDIVB *eta_KDIVB *Rc_ENG_KDIVB;   ,
            scrh += alambda_KSLOWM*eta_KSLOWM*Rc_ENG_KSLOWM; 
            scrh += alambda_KSLOWP*eta_KSLOWP*Rc_ENG_KSLOWP;  ,
            scrh += alambda_KALFVM*eta_KALFVM*Rc_ENG_KALFVM; 
            scrh += alambda_KALFVP*eta_KALFVP*Rc_ENG_KALFVP;  )
        scrh += alambda_KENTRP*eta_KENTRP*Rc_ENG_KENTRP; 
        Flux(ENG,k,j,i) = HALF_F*(fluxL_ENG + fluxR_ENG - scrh);
#endif
            
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


    Kokkos::Profiling::popRegion();

}
