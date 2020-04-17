#include "../idefix.hpp"
#include "solvers.hpp"

// Compute Riemann fluxes from states using TVDLF solver
void Tvdlf(DataBlock & data, int dir, real gamma, real C2Iso) {
    int ioffset,joffset,koffset;

    Kokkos::Profiling::pushRegion("TVDLF_Solver");
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

    real gamma_m1=gamma-ONE_F;

    idefix_for("TVDLF_Kernel",data.beg[KDIR],data.end[KDIR]+koffset,data.beg[JDIR],data.end[JDIR]+joffset,data.beg[IDIR],data.end[IDIR]+ioffset,
                        KOKKOS_LAMBDA (int k, int j, int i) 
            {
                int VXn = VX1+dir;
                // Primitive variables
                real vL[NVAR];
                real vR[NVAR];
                real vVXn;

                // temporary array to compute conservative variables and fluxes
                real u[NVAR];

                // temporary flux array
                real fl[NVAR];

                // Signal speeds
                real cRL, cmax;

                // 1-- Store the primitive variables on the left, right, and averaged states
                for(int nv = 0 ; nv < NVAR; nv++) {
                    vL[nv] = PrimL(nv,k,j,i);
                    vR[nv] = PrimR(nv,k,j,i);
                }
                vVXn = HALF_F*(vL[VXn]+vR[VXn]);
            #if HAVE_ENERGY
                real vRHO = HALF_F*(vL[RHO]+vR[RHO]);
                real vPRS = HALF_F*(vL[PRS]+vR[PRS]);
            #endif
                
                // 2-- Get the wave speed
            #if HAVE_ENERGY
                real a2 = vPRS / vRHO;
                cRL = SQRT(a2*(gamma_m1+ONE_F));
            #else
                cRL = SQRT(C2Iso);
            #endif
                cmax = FMAX(FABS(vVXn+cRL),FABS(vVXn-cRL));

                // 3-- Compute the conservative variables on the left
                K_PrimToCons(u, vL, gamma_m1);
                // store result in flux array
                for(int nv = 0 ; nv < NVAR; nv++) {
                    fl[nv] = HALF_F*cmax*u[nv];
                }

                // 4-- Compute the left fluxes
                K_Flux(u, vL, u, C2Iso, dir);

                // 5-- Compute the flux from the left state
                for(int nv = 0 ; nv < NVAR; nv++) {
                    fl[nv] += HALF_F*u[nv];
                }

                // 6-- Compute the conservative variables on the right
                K_PrimToCons(u, vR, gamma_m1);
                // store uL
                for(int nv = 0 ; nv < NVAR; nv++) {
                    fl[nv] -= HALF_F*cmax*u[nv];
                }

                // 7-- Compute the right fluxes
                K_Flux(u, vR, u, C2Iso, dir);

                // 8-- Compute the flux from the left and right states
                for(int nv = 0 ; nv < NVAR; nv++) {
                    Flux(nv,k,j,i) = fl[nv] + HALF_F*u[nv];
                }

                //9-- Compute maximum dt for this sweep
                const int ig = ioffset*i + joffset*j + koffset*k;

                invDt(k,j,i) += cmax/dx(ig);

            });

    Kokkos::Profiling::popRegion();

}

