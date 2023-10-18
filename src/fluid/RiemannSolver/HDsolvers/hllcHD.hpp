// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_RIEMANNSOLVER_HDSOLVERS_HLLCHD_HPP_
#define FLUID_RIEMANNSOLVER_HDSOLVERS_HLLCHD_HPP_

#include "../idefix.hpp"
#include "fluid.hpp"
#include "extrapolateToFaces.hpp"
#include "flux.hpp"
#include "convertConsToPrim.hpp"

// Compute Riemann fluxes from states using HLLC solver
template <typename Phys>
template<const int DIR>
void RiemannSolver<Phys>::HllcHD(IdefixArray4D<real> &Flux) {
  idfx::pushRegion("RiemannSolver::HLLC_Solver");

  constexpr int ioffset = (DIR==IDIR) ? 1 : 0;
  constexpr int joffset = (DIR==JDIR) ? 1 : 0;
  constexpr int koffset = (DIR==KDIR) ? 1 : 0;

  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Vs = this->Vs;
  IdefixArray3D<real> cMax = this->cMax;

  EquationOfState eos = *(hydro->eos.get());

  ExtrapolateToFaces<Phys,DIR> extrapol = *this->GetExtrapolator<DIR>();

  idefix_for("HLLC_Kernel",
             data->beg[KDIR],data->end[KDIR]+koffset,
             data->beg[JDIR],data->end[JDIR]+joffset,
             data->beg[IDIR],data->end[IDIR]+ioffset,
    KOKKOS_LAMBDA (int k, int j, int i) {
      // Init the directions (should be in the kernel for proper optimisation by the compilers)
      EXPAND( constexpr int Xn = DIR+MX1;                    ,
              constexpr int Xt = (DIR == IDIR ? MX2 : MX1);  ,
              constexpr int Xb = (DIR == KDIR ? MX2 : MX3);  )

      // Primitive variables
      real vL[Phys::nvar];
      real vR[Phys::nvar];

      // Conservative variables
      real uL[Phys::nvar];
      real uR[Phys::nvar];

      // Flux (left and right)
      real fluxL[Phys::nvar];
      real fluxR[Phys::nvar];

      // Signal speeds
      real cL, cR, cmax;

      // 1-- Store the primitive variables on the left, right, and averaged states
      extrapol.ExtrapolatePrimVar(i, j, k, vL, vR);

      // 2-- Get the wave speed
      #if HAVE_ENERGY
        cL = std::sqrt(eos.GetGamma(vL[PRS],vL[RHO])*(vL[PRS]/vL[RHO]));
        cR = std::sqrt(eos.GetGamma(vR[PRS],vR[RHO])*(vR[PRS]/vR[RHO]));
      #else
        cL = HALF_F*(eos.GetWaveSpeed(k,j,i)
                    +eos.GetWaveSpeed(k-koffset,j-joffset,i-ioffset));
        cR = cL;
      #endif

      real cminL = vL[Xn] - cL;
      real cmaxL = vL[Xn] + cL;

      real cminR = vR[Xn] - cR;
      real cmaxR = vR[Xn] + cR;

      real SL = FMIN(cminL, cminR);
      real SR = FMAX(cmaxL, cmaxR);

      cmax  = FMAX(FABS(SL), FABS(SR));

      // 3-- Compute the conservative variables
      K_PrimToCons<Phys>(uL, vL, &eos);
      K_PrimToCons<Phys>(uR, vR, &eos);

      // 4-- Compute the left and right fluxes
      K_Flux<Phys,DIR>(fluxL, vL, uL, cL*cL);
      K_Flux<Phys,DIR>(fluxR, vR, uR, cR*cR);

      // 5-- Compute the flux from the left and right states
      if (SL > 0) {
#pragma unroll
        for (int nv = 0 ; nv < Phys::nvar; nv++) {
          Flux(nv,k,j,i) = fluxL[nv];
        }
      } else if (SR < 0) {
#pragma unroll
        for (int nv = 0 ; nv < Phys::nvar; nv++) {
          Flux(nv,k,j,i) = fluxR[nv];
        }
      } else {
        real usL[Phys::nvar];
        real usR[Phys::nvar];
        real vs;

#if HAVE_ENERGY
        real qL, qR, wL, wR;
        qL = vL[PRS] + uL[Xn]*(vL[Xn] - SL);
        qR = vR[PRS] + uR[Xn]*(vR[Xn] - SR);

        wL = vL[RHO]*(vL[Xn] - SL);
        wR = vR[RHO]*(vR[Xn] - SR);

        vs = (qR - qL)/(wR - wL); // wR - wL > 0 since SL < 0, SR > 0

        usL[RHO] = uL[RHO]*(SL - vL[Xn])/(SL - vs);
        usR[RHO] = uR[RHO]*(SR - vR[Xn])/(SR - vs);
        EXPAND(usL[Xn] = usL[RHO]*vs;     usR[Xn] = usR[RHO]*vs;      ,
                usL[Xt] = usL[RHO]*vL[Xt]; usR[Xt] = usR[RHO]*vR[Xt];  ,
                usL[Xb] = usL[RHO]*vL[Xb]; usR[Xb] = usR[RHO]*vR[Xb];)

        usL[ENG] =    uL[ENG]/vL[RHO]
                    + (vs - vL[Xn])*(vs + vL[PRS]/(vL[RHO]*(SL - vL[Xn])));
        usR[ENG] =    uR[ENG]/vR[RHO]
                    + (vs - vR[Xn])*(vs + vR[PRS]/(vR[RHO]*(SR - vR[Xn])));

        usL[ENG] *= usL[RHO];
        usR[ENG] *= usR[RHO];
#else
        real scrh = 1.0/(SR - SL);
        real rho  = (SR*uR[RHO] - SL*uL[RHO] - fluxR[RHO] + fluxL[RHO])*scrh;
        real mx   = (SR*uR[Xn] - SL*uL[Xn] - fluxR[Xn] + fluxL[Xn])*scrh;

        usL[RHO] = usR[RHO] = rho;
        usL[Xn] = usR[Xn] = mx;
        vs  = (  SR*fluxL[RHO] - SL*fluxR[RHO]
                + SR*SL*(uR[RHO] - uL[RHO]));
        vs *= scrh;
        vs /= rho;
        EXPAND(                                            ,
                usL[Xt] = rho*vL[Xt]; usR[Xt] = rho*vR[Xt]; ,
                usL[Xb] = rho*vL[Xb]; usR[Xb] = rho*vR[Xb];)
#endif

    // Compute the flux from the left and right states
        if (vs >= 0.0) {
#pragma unroll
          for(int nv = 0 ; nv < Phys::nvar; nv++) {
            Flux(nv,k,j,i) = fluxL[nv] + SL*(usL[nv] - uL[nv]);
          }
        } else {
#pragma unroll
          for(int nv = 0 ; nv < Phys::nvar; nv++) {
            Flux(nv,k,j,i) = fluxR[nv] + SR*(usR[nv] - uR[nv]);
          }
        }
      }

      //6-- Compute maximum wave speed for this sweep
      cMax(k,j,i) = cmax;
  });

  idfx::popRegion();
}

#endif  // FLUID_RIEMANNSOLVER_HDSOLVERS_HLLCHD_HPP_
