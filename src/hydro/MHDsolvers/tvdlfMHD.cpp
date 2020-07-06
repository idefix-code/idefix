#include "../idefix.hpp"
#include "solversMHD.hpp"

// Compute Riemann fluxes from states using TVDLF solver
void TvdlfMHD(DataBlock & data, int dir, real gamma, real C2Iso) {
    int ioffset,joffset,koffset;
    int iextend, jextend,kextend;

    idfx::pushRegion("TVDLF_MHD");
    
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
        
        real v_RHO;
        EXPAND(
            real v_VX1;
            real v_BX1; ,
            real v_VX2;
            real v_BX2; ,
            real v_VX3;
            real v_BX3; )

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

        // 1.3-- compute primitive variables in the center
        v_RHO = HALF_F*(vL_RHO + vR_RHO);
        EXPAND(
            v_VX1 = HALF_F*(vL_VX1 + vR_VX1);
            v_BX1 = HALF_F*(vL_BX1 + vR_BX1); ,
            v_VX2 = HALF_F*(vL_VX2 + vR_VX2);
            v_BX2 = HALF_F*(vL_BX2 + vR_BX2); ,
            v_VX3 = HALF_F*(vL_VX3 + vR_VX3);
            v_BX3 = HALF_F*(vL_BX3 + vR_BX3); )
        
        EXPAND(
            real v_VXn;
            real v_BXn; ,
            real v_BXt; ,
            real v_BXb; )
        
        EXPAND (
        if (dir == IDIR) {
            EXPAND (
            v_VXn = v_VX1;
            v_BXn = v_BX1;  ,
            v_BXt = v_BX2;  ,
            v_BXb = v_BX3;  )
        }
        ,
        if (dir == JDIR) {
            EXPAND (
            v_VXn = v_VX2;
            v_BXn = v_BX2;  ,
            v_BXt = v_BX1;  ,
            v_BXb = v_BX3;  )
        }
        ,
        if (dir == KDIR) {
            EXPAND (
            v_VXn = v_VX3;
            v_BXn = v_BX3;  ,
            v_BXt = v_BX1;  ,
            v_BXb = v_BX2;  )
        }
        )

#if HAVE_ENERGY
            real vL_PRS = PrimL(PRS,k,j,i);
            real vR_PRS = PrimR(PRS,k,j,i);
            real v_PRS = HALF_F*(vL_PRS + vR_PRS);
            real uL_ENG, uR_ENG;
            real fluxL_ENG, fluxR_ENG;
#endif
            
            
        // 2-- Get the wave speed
        real c, cmax;
        real gpr, Bt2, B2;

        #if HAVE_ENERGY
            gpr=gamma*v_PRS;
        #else
            gpr=C2Iso*v_RHO;
        #endif
        Bt2=EXPAND(ZERO_F    ,
                    + v_BXt*v_BXt,
                    + v_BXb*v_BXb);

        B2=Bt2 + v_BXn*v_BXn;

        c = gpr - B2;
        c = c + B2 + SQRT(c*c + FOUR_F*gpr*Bt2);
        c = SQRT(HALF_F * c/v_RHO);

        cmax = FMAX(FABS(v_VXn+c),FABS(v_VXn-c));
        
        
        // 2-- Compute the conservative variables
        //K_PrimToCons(u, v, gamma_m1);
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
        Flux(RHO,k,j,i) = HALF_F*(fluxL_RHO+fluxR_RHO - cmax*(uR_RHO-uL_RHO));
        EXPAND ( Flux(MX1,k,j,i) = HALF_F * (fluxL_MX1 + fluxR_MX1 - cmax*(uR_MX1-uL_MX1));
                 Flux(BX1,k,j,i) = HALF_F * (fluxL_BX1 + fluxR_BX1 - cmax*(uR_BX1-uL_BX1));    ,
                 Flux(MX2,k,j,i) = HALF_F * (fluxL_MX2 + fluxR_MX2 - cmax*(uR_MX2-uL_MX2));
                 Flux(BX2,k,j,i) = HALF_F * (fluxL_BX2 + fluxR_BX2 - cmax*(uR_BX2-uL_BX2));    ,
                 Flux(MX3,k,j,i) = HALF_F * (fluxL_MX3 + fluxR_MX3 - cmax*(uR_MX3-uL_MX3));
                 Flux(BX3,k,j,i) = HALF_F * (fluxL_BX3 + fluxR_BX3 - cmax*(uR_BX3-uL_BX3));    )
#if HAVE_ENERGY
        Flux(ENG,k,j,i) = HALF_F*(fluxL_ENG+fluxR_ENG - cmax*(uR_ENG-uL_ENG));
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


    idfx::popRegion();

}

