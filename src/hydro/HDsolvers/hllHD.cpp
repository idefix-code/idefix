#include "../idefix.hpp"
#include "solversHD.hpp"

// Compute Riemann fluxes from states using HLL solver
void HllHD(DataBlock & data, int dir, real gamma, real C2Iso) {
    int ioffset,joffset,koffset;

    Kokkos::Profiling::pushRegion("HLL_Solver");
    ioffset=joffset=koffset=0;
    // Determine the offset along which we do the extrapolation
    if(dir==IDIR) ioffset=1;
    if(dir==JDIR) joffset=1;
    if(dir==KDIR) koffset=1;

    IdefixArray4D<real> PrimL = data.PrimL;
    IdefixArray4D<real> PrimR = data.PrimR;
    IdefixArray4D<real> Flux = data.FluxRiemann;
    IdefixArray1D<real> dx = data.dx[dir];
    IdefixArray3D<real> invDt = data.InvDtHyp;

    real gamma_m1 = gamma - ONE_F;

    idefix_for("HLL_Kernel",data.beg[KDIR],data.end[KDIR]+koffset,data.beg[JDIR],data.end[JDIR]+joffset,data.beg[IDIR],data.end[IDIR]+ioffset,
                        KOKKOS_LAMBDA (int k, int j, int i) 
            {
                int VXn = VX1+dir;

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
#if HAVE_ENERGY
                cL = SQRT( gamma *(vL[PRS]/vL[RHO]));
                cR = SQRT( gamma *(vR[PRS]/vR[RHO]));
#else
                cL = SQRT(C2Iso);
                cR = cL;
#endif
                
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
                K_Flux(fluxL, vL, fluxL, C2Iso, dir);
                K_Flux(fluxR, vR, fluxR, C2Iso, dir);

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

                invDt(k,j,i) += cmax/dx(ig);

            });

    Kokkos::Profiling::popRegion();

}
