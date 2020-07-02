#include "../idefix.hpp"
#include "solversMHD.hpp"

// Compute Riemann fluxes from states using HLL solver
void HllMHD(DataBlock & data, int dir, real gamma, real C2Iso) {

    int ioffset,joffset,koffset;
    int iextend, jextend,kextend;

    Kokkos::Profiling::pushRegion("HLL_MHD");
    
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
    // st and sb will be useful only when Hall is included
    EXPAND( int BXt;  ,
                      ,
            int BXb;  )
    
    D_EXPAND( real st;  ,
                        ,
              real sb;  )

    switch(dir) {
        case(IDIR):
            ioffset = 1;
            D_EXPAND(               ,
                   jextend = 1;     ,
                   kextend = 1; )

            EXPAND(                 , 
                   BXt = BX2;       , 
                   BXb = BX3;       )

            Et = data.emf.ezi;
            Eb = data.emf.eyi;
            SV = data.emf.svx;

            D_EXPAND( st = -1.0;  ,
                                  ,
                      sb = +1.0;  )
            break;
        case(JDIR):
            joffset=1;
            D_EXPAND( iextend = 1;  ,
                                    ,
                    kextend = 1;)
            EXPAND(                 , 
                   BXt = BX1;       , 
                   BXb = BX3;       )

            Et = data.emf.ezj;
            Eb = data.emf.exj;
            SV = data.emf.svy;

            D_EXPAND( st = +1.0;  ,
                                  ,
                      sb = -1.0;  )
            break;
        case(KDIR):
            koffset=1;
            D_EXPAND( iextend = 1;  ,
                    jextend = 1;    ,
                    )
            EXPAND(                 , 
                   BXt = BX1;       , 
                   BXb = BX2;       )

            Et = data.emf.eyk;
            Eb = data.emf.exk;
            SV = data.emf.svz;

            D_EXPAND( st = -1.0;  ,
                                  ,
                      sb = +1.0;  )
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
            Flux(RHO,k,j,i) = (SL*SR*uR_RHO - SL*SR*uL_RHO + SR*fluxL_RHO - SL*fluxR_RHO) / (SR - SL);
            EXPAND ( Flux(MX1,k,j,i) = (SL*SR*uR_MX1 - SL*SR*uL_MX1 + SR*fluxL_MX1 - SL*fluxR_MX1) / (SR - SL);
                 Flux(BX1,k,j,i) = (SL*SR*uR_BX1 - SL*SR*uL_BX1 + SR*fluxL_BX1 - SL*fluxR_BX1) / (SR - SL);    ,
                 Flux(MX2,k,j,i) = (SL*SR*uR_MX2 - SL*SR*uL_MX2 + SR*fluxL_MX2 - SL*fluxR_MX2) / (SR - SL);
                 Flux(BX2,k,j,i) = (SL*SR*uR_BX2 - SL*SR*uL_BX2 + SR*fluxL_BX2 - SL*fluxR_BX2) / (SR - SL);    ,
                 Flux(MX3,k,j,i) = (SL*SR*uR_MX3 - SL*SR*uL_MX3 + SR*fluxL_MX3 - SL*fluxR_MX3) / (SR - SL);
                 Flux(BX3,k,j,i) = (SL*SR*uR_BX3 - SL*SR*uL_BX3 + SR*fluxL_BX3 - SL*fluxR_BX3) / (SR - SL);    )
#if HAVE_ENERGY
            Flux(ENG,k,j,i) = (SL*SR*uR_ENG - SL*SR*uL_ENG + SR*fluxL_ENG - SL*fluxR_ENG) / (SR - SL);
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
        if      (Flux(RHO,k,j,i) >  eps_UCT_CONTACT) s =  1;
        else if (Flux(RHO,k,j,i) < -eps_UCT_CONTACT) s = -1;

        SV(k,j,i) = s;
#endif

    });


    Kokkos::Profiling::popRegion();

}
