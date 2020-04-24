#include "../idefix.hpp"
#include "solvers.hpp"

// Compute Riemann fluxes from states using TVDLF solver
void Tvdlf(DataBlock & data, int dir, real gamma, real C2Iso) {
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
        real v[NVAR];
        real u[NVAR];
        real flux[NVAR];
        real fluxRiemann[NVAR];

        // Store the average primitive variables
        for(int nv = 0 ; nv < NVAR; nv++) {
            v[nv] = HALF_F*(PrimL(nv,k,j,i) + PrimR(nv,k,j,i));
        }


        // 4-- Get the wave speed
        // Signal speeds
        real cRL, cmax;
        real gpr, Bt2, B2;

        #if HAVE_ENERGY
            gpr=(gamma_m1+ONE_F)*v[PRS];
        #else
            gpr=C2Iso*v[RHO];
        #endif
        Bt2=EXPAND(ZERO_F    ,
                    + v[BXt]*v[BXt],
                    + v[BXb]*v[BXb]);

        B2=Bt2 + v[BXn]*v[BXn];

        cRL = gpr - B2;
        cRL = cRL + B2 + SQRT(cRL*cRL + FOUR_F*gpr*Bt2);
        cRL = SQRT(HALF_F * cRL/v[RHO]);

        cmax = FMAX(FABS(v[VXn]+cRL),FABS(v[VXn]-cRL));


        // Load the left state
        for(int nv = 0 ; nv < NVAR; nv++) {
            v[nv] = PrimL(nv,k,j,i);
        }

        // 2-- Compute the conservative variables
        K_PrimToCons(u, v, gamma_m1);

        // 3-- Compute the left and right fluxes
        K_Flux(flux, v, u, C2Iso, VXn, VXt, VXb, BXn, BXt, BXb, MXn);
        

        // 5-- Compute the flux from the left and right states
        for(int nv = 0 ; nv < NVAR; nv++) {
            fluxRiemann[nv] = flux[nv] + cmax*u[nv];
        }

        // Load the right state
        for(int nv = 0 ; nv < NVAR; nv++) {
            v[nv] = PrimR(nv,k,j,i);
        }

        // 2-- Compute the conservative variables
        K_PrimToCons(u, v, gamma_m1);

        // 3-- Compute the left and right fluxes
        K_Flux(flux, v, u, C2Iso, VXn, VXt, VXb, BXn, BXt, BXb, MXn);
        
        // 5-- Compute the flux from the left and right states
        for(int nv = 0 ; nv < NVAR; nv++) {
            Flux(nv,k,j,i) = HALF_F*(fluxRiemann[nv]+flux[nv] - cmax*u[nv]);
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

