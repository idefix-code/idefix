#include "../idefix.hpp"

// Compute Riemann fluxes from states using HLL solver
template<const int DIR, const int Xn, const int Xt, const int Xb>
void HllHD(DataBlock & data, real gamma, real C2Iso) {
    int ioffset,joffset,koffset;

    idfx::pushRegion("HLL_Solver");
    ioffset=joffset=koffset=0;
    // Determine the offset along which we do the extrapolation
    if(DIR==IDIR) ioffset=1;
    if(DIR==JDIR) joffset=1;
    if(DIR==KDIR) koffset=1;

    IdefixArray4D<real> PrimL = data.PrimL;
    IdefixArray4D<real> PrimR = data.PrimR;
    IdefixArray4D<real> Flux = data.FluxRiemann;
    IdefixArray3D<real> cMax = data.cMax;

    real gamma_m1 = gamma - ONE_F;

    idefix_for("HLL_Kernel",data.beg[KDIR],data.end[KDIR]+koffset,data.beg[JDIR],data.end[JDIR]+joffset,data.beg[IDIR],data.end[IDIR]+ioffset,
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
                #pragma unroll
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
                real cminL = vL[Xn] - cL;
                real cmaxL = vL[Xn] + cL;
                
                real cminR = vR[Xn] - cR;
                real cmaxR = vR[Xn] + cR;
                
                real SL = FMIN(cminL, cminR);
                real SR = FMAX(cmaxL, cmaxR);
                
                cmax  = FMAX(FABS(SL), FABS(SR));
                
                // 2-- Compute the conservative variables
                K_PrimToCons(uL, vL, gamma_m1);
                K_PrimToCons(uR, vR, gamma_m1);

                // 3-- Compute the left and right fluxes
                K_Flux(fluxL, vL, uL, C2Iso, Xn);
                K_Flux(fluxR, vR, uR, C2Iso, Xn);

                // 5-- Compute the flux from the left and right states
                if (SL > 0){
                    #pragma unroll
                    for (int nv = 0 ; nv < NFLX; nv++) {
                        Flux(nv,k,j,i) = fluxL[nv];
                    }
                }
                else if (SR < 0) {
                    #pragma unroll
                    for (int nv = 0 ; nv < NFLX; nv++) {
                        Flux(nv,k,j,i) = fluxR[nv];
                    }
                }
                else {
                    #pragma unroll
                    for(int nv = 0 ; nv < NFLX; nv++) {
                        Flux(nv,k,j,i) = SL*SR*uR[nv] - SL*SR*uL[nv] + SR*fluxL[nv] - SL*fluxR[nv];
                        Flux(nv,k,j,i) /= (SR - SL);
                    }
                }

                //6-- Compute maximum wave speed for this sweep

                cMax(k,j,i) = cmax;

            });

    idfx::popRegion();

}
