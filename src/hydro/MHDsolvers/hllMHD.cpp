#include "../idefix.hpp"
#include "solversMHD.hpp"

// Compute Riemann fluxes from states using HLL solver
void HllMHD(DataBlock & data, int dir, real gamma, real C2Iso) {

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
        
        for(int nv = 0 ; nv < NVAR; nv++) {
            fluxL[nv] = uL[nv];
            fluxR[nv] = uR[nv];
        }

        // 3-- Compute the left and right fluxes
        K_Flux(fluxL, vL, fluxL, C2Iso, VXn, VXt, VXb, BXn, BXt, BXb, MXn);
        K_Flux(fluxR, vR, fluxR, C2Iso, VXn, VXt, VXb, BXn, BXt, BXb, MXn);

        // 5-- Compute the flux from the left and right states
        if (SL > 0){
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
            for(int nv = 0 ; nv < NFLX; nv++) {
                Flux(nv,k,j,i) = SL*SR*uR[nv] - SL*SR*uL[nv] + SR*fluxL[nv] - SL*fluxR[nv];
                Flux(nv,k,j,i) *= (1.0 / (SR - SL));
            }
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
