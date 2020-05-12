#include "../idefix.hpp"
#include "../solversHD.hpp"

// Compute Riemann fluxes from states using HLLC solver
void HllcHD(DataBlock & data, int dir, real gamma, real C2Iso) {
    int ioffset,joffset,koffset;

    Kokkos::Profiling::pushRegion("HLLC_Solver");
    ioffset=joffset=koffset=0;
    // Determine the offset along which we do the extrapolation
    EXPAND( int VXn; int MXn;  ,
            int VXt; int MXt;  ,
            int VXb; int MXb;   )
    
    switch(dir) {
        case(IDIR):
            ioffset = 1;
            
            EXPAND(VXn = MXn = VX1;  , 
                   VXt = MXt = VX2;  , 
                   VXb = MXb = VX3; )
            break;
        case(JDIR):
            joffset=1;

            EXPAND(VXn = MXn = VX2;  , 
                   VXt = MXt = VX1;  , 
                   VXb = MXb = VX3; )
            break;
        case(KDIR):
            koffset=1;

            EXPAND(VXn = MXn = VX3;  , 
                   VXt = MXt = VX1;  , 
                   VXb = MXb = VX2; )
            break;
        default:
            IDEFIX_ERROR("Wrong direction");
    }

    IdefixArray4D<real> PrimL = data.PrimL;
    IdefixArray4D<real> PrimR = data.PrimR;
    IdefixArray4D<real> Flux = data.FluxRiemann;
    IdefixArray1D<real> dx = data.dx[dir];
    IdefixArray3D<real> invDt = data.InvDtHyp;

    real gamma_m1 = gamma - ONE_F;

    idefix_for("HLLC_Kernel",data.beg[KDIR],data.end[KDIR]+koffset,data.beg[JDIR],data.end[JDIR]+joffset,data.beg[IDIR],data.end[IDIR]+ioffset,
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

                // 3-- Compute the left and right fluxes
                for(int nv = 0 ; nv < NVAR; nv++) {
                    fluxL[nv] = uL[nv];
                    fluxR[nv] = uR[nv];
                }
                
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
                    real usL[NVAR];
                    real usR[NVAR];
                    real vs;
                    
#if HAVE_ENERGY
                    real qL, qR, wL, wR;
                    qL = vL[PRS] + uL[MXn]*(vL[VXn] - SL);
                    qR = vR[PRS] + uR[MXn]*(vR[VXn] - SR);

                    wL = vL[RHO]*(vL[VXn] - SL);
                    wR = vR[RHO]*(vR[VXn] - SR);

                    vs = (qR - qL)/(wR - wL); // wR - wL > 0 since SL < 0, SR > 0

                    usL[RHO] = uL[RHO]*(SL - vL[VXn])/(SL - vs);
                    usR[RHO] = uR[RHO]*(SR - vR[VXn])/(SR - vs);
                    EXPAND(usL[MXn] = usL[RHO]*vs;     usR[MXn] = usR[RHO]*vs;      ,
                           usL[MXt] = usL[RHO]*vL[VXt]; usR[MXt] = usR[RHO]*vR[VXt];  ,
                           usL[MXb] = usL[RHO]*vL[VXb]; usR[MXb] = usR[RHO]*vR[VXb];)

                    usL[ENG] =    uL[ENG]/vL[RHO]
                               + (vs - vL[VXn])*(vs + vL[PRS]/(vL[RHO]*(SL - vL[VXn])));
                    usR[ENG] =    uR[ENG]/vR[RHO]
                               + (vs - vR[VXn])*(vs + vR[PRS]/(vR[RHO]*(SR - vR[VXn])));

                    usL[ENG] *= usL[RHO];
                    usR[ENG] *= usR[RHO];
#else
                    real scrh = 1.0/(SR - SL);
                    real rho  = (SR*uR[RHO] - SL*uL[RHO] - fluxR[RHO] + fluxL[RHO])*scrh;
                    real mx   = (SR*uR[MXn] - SL*uL[MXn] - fluxR[MXn] + fluxL[MXn])*scrh;

                    usL[RHO] = usR[RHO] = rho;
                    usL[MXn] = usR[MXn] = mx;
                    vs  = (  SR*fluxL[RHO] - SL*fluxR[RHO]
                           + SR*SL*(uR[RHO] - uL[RHO]));
                    vs *= scrh;
                    vs /= rho;
                    EXPAND(                                            ,
                           usL[MXt] = rho*vL[VXt]; usR[MXt] = rho*vR[VXt]; ,
                           usL[MXb] = rho*vL[VXb]; usR[MXb] = rho*vR[VXb];)
#endif
                    
                // 5-- Compute the flux from the left and right states
                    if (vs >= 0.0){
                        for(int nv = 0 ; nv < NVAR; nv++) {
                            Flux(nv,k,j,i) = fluxL[nv] + SL*(usL[nv] - uL[nv]);
                        }
                    } else {
                        for(int nv = 0 ; nv < NVAR; nv++) {
                            Flux(nv,k,j,i) = fluxR[nv] + SR*(usR[nv] - uR[nv]);
                        }
                    }
                }

                //6-- Compute maximum dt for this sweep
                const int ig = ioffset*i + joffset*j + koffset*k;

                invDt(k,j,i) += cmax/dx(ig);
            });

    Kokkos::Profiling::popRegion();

}
