// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_RIEMANNSOLVER_DUSTSOLVERS_HLLDUST_HPP_
#define FLUID_RIEMANNSOLVER_DUSTSOLVERS_HLLDUST_HPP_

#include "../idefix.hpp"
#include "fluid.hpp"
#include "slopeLimiter.hpp"
#include "flux.hpp"
#include "convertConsToPrim.hpp"

// Compute Riemann fluxes from states using HLL solver
template <typename Phys>
template<const int DIR>
void RiemannSolver<Phys>::HllDust(IdefixArray4D<real> &Flux) {
  idfx::pushRegion("Hydro::HLL_Dust");

  constexpr int ioffset = (DIR==IDIR) ? 1 : 0;
  constexpr int joffset = (DIR==JDIR) ? 1 : 0;
  constexpr int koffset = (DIR==KDIR) ? 1 : 0;

  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray3D<real> cMax = this->cMax;

  // Required for high order interpolations
  IdefixArray1D<real> dx = this->data->dx[DIR];


  SlopeLimiter<Phys,DIR> slopeLim(Vc,data->dx[DIR],haveShockFlattening,shockFlattening.get());;

  idefix_for("HLL_Kernel",
             data->beg[KDIR],data->end[KDIR]+koffset,
             data->beg[JDIR],data->end[JDIR]+joffset,
             data->beg[IDIR],data->end[IDIR]+ioffset,
    KOKKOS_LAMBDA (int k, int j, int i) {
      // Init the directions (should be in the kernel for proper optimisation by the compilers)
      constexpr int Xn = DIR+MX1;

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
      slopeLim.ExtrapolatePrimVar(i, j, k, vL, vR);

      // 2-- Get the wave speed

      real SL = vL[Xn];
      real SR = vR[Xn];

      cmax  = FMAX(FABS(SL), FABS(SR));

      // 3-- Compute the conservative variables: do this by extrapolation
      K_PrimToCons<Phys>(uL, vL, 0.0); // Set gamma to 0 implicitly
      K_PrimToCons<Phys>(uR, vR, 0.0);

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
        real dS = SR-SL;
        if(std::abs(dS) < SMALL_NUMBER) {
          dS = SMALL_NUMBER;
        }
#pragma unroll
        for(int nv = 0 ; nv < Phys::nvar; nv++) {
          Flux(nv,k,j,i) = SL*SR*uR[nv] - SL*SR*uL[nv] + SR*fluxL[nv] - SL*fluxR[nv];
          Flux(nv,k,j,i) /= dS;
        }
      }

      //6-- Compute maximum wave speed for this sweep
      cMax(k,j,i) = cmax;
    }
  );

  idfx::popRegion();
}

#endif // FLUID_RIEMANNSOLVER_DUSTSOLVERS_HLLDUST_HPP_
