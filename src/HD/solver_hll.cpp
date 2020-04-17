#include "../idefix.hpp"
#include "solvers.hpp"

// Compute Riemann fluxes from states using HLL solver
void Hll(DataBlock & data, int dir, real gamma, real C2Iso) {
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
                int MXn = VXn;
                // Primitive variables
                real vL[NVAR];
                real vR[NVAR];
                real vRL[NVAR];

                // Conservative variables
                real uL[NVAR];
                real uR[NVAR];

                // Flux (left and right)
                real fluxL[NVAR];
                real fluxR[NVAR];

                // Signal speeds
                real cRL, cmax;
                real SR, SL, sl_min, sr_min, sl_max, sr_max;
                real aL, a2L, a2R, aR;
                real scrh;

                // 1-- Store the primitive variables on the left, right, and averaged states
                for(int nv = 0 ; nv < NVAR; nv++) {
                    vL[nv] = PrimL(nv,k,j,i);
                    vR[nv] = PrimR(nv,k,j,i);
                    vRL[nv]=HALF_F*(vL[nv]+vR[nv]);
                }

                // 2-- Compute the conservative variables
                K_PrimToCons(uL, vL, gamma_m1);
                K_PrimToCons(uR, vR, gamma_m1);

                // 3-- Compute the left and right fluxes
                K_Flux(fluxL, vL, uL, C2Iso, dir);
                K_Flux(fluxR, vR, uR, C2Iso, dir);

                // 4-- Get the wave speed
                
                // HLL_Speed using DAVIS_ESTIMATE
                a2L = gamma*vL[PRS]/vL[RHO];
                a2R = gamma*vR[PRS]/vR[RHO];
                
                aL = sqrt(a2L);
                aR = sqrt(a2R);

                sl_min = vL[VXn] - aL;
                sl_max = vL[VXn] + aL;
                
                sr_min = vR[VXn] - aR;
                sr_max = vR[VXn] + aR;

                SL = FMIN(sl_min, sr_min);
                SR = FMAX(sl_max, sr_max);

                /*scrh  = fabs(vL[VXn]) + fabs(vR[VXn]);    
                scrh /= aL + aR;

                // reduction!!!
                g_maxMach = MAX(scrh, g_maxMach);*/
                
                cmax  = FMAX(fabs(SL), fabs(SR));

                // 5-- Compute the flux from the left and right states
                if (SL > 0.0) {
                    for(int nv = 0 ; nv < NVAR; nv++) Flux(nv,k,j,i) = fluxL[nv];
                }
                else if (SR < 0.0) {
                    for(int nv = 0 ; nv < NVAR; nv++) Flux(nv,k,j,i) = fluxR[nv];
                }
                else {
                    scrh = 1.0 / (SR - SL);
                    for(int nv = 0 ; nv < NVAR; nv++) {
                        Flux(nv,k,j,i) = SL*SR*(uR[nv] - uL[nv]) +
                                         SR*fluxL[nv] - SL*fluxR[nv];
                        Flux(nv,k,j,i) *= scrh;
                    }
                }

                //6-- Compute maximum dt for this sweep
                const int ig = ioffset*i + joffset*j + koffset*k;
                invDt(k,j,i) += cmax/dx(ig);
            });

    Kokkos::Profiling::popRegion();

}
