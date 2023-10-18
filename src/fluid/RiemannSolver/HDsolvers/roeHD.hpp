// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_RIEMANNSOLVER_HDSOLVERS_ROEHD_HPP_
#define FLUID_RIEMANNSOLVER_HDSOLVERS_ROEHD_HPP_

#include "../idefix.hpp"
#include "fluid.hpp"
#include "extrapolateToFaces.hpp"
#include "flux.hpp"
#include "convertConsToPrim.hpp"

#define ROE_AVERAGE 0
#undef NMODES

#if HAVE_ENERGY
  #define I0 0
  #define I1 1
  #define IE 2
  #define I2 3
  #define I3 4

  #define NMODES 5

#else
  #define I0 0
  #define I1 1
  #define I2 2
  #define I3 3

  #define NMODES 4
#endif

// Compute Riemann fluxes from states using ROE solver
template <typename Phys>
template<const int DIR>
void RiemannSolver<Phys>::RoeHD(IdefixArray4D<real> &Flux) {
  idfx::pushRegion("RiemannSolver::ROE_Solver");

  constexpr int ioffset = (DIR==IDIR) ? 1 : 0;
  constexpr int joffset = (DIR==JDIR) ? 1 : 0;
  constexpr int koffset = (DIR==KDIR) ? 1 : 0;

  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Vs = this->Vs;
  IdefixArray3D<real> cMax = this->cMax;

  EquationOfState eos = *(hydro->eos.get());

  ExtrapolateToFaces<Phys,DIR> extrapol = *this->GetExtrapolator<DIR>();

  idefix_for("ROE_Kernel",
             data->beg[KDIR],data->end[KDIR]+koffset,
             data->beg[JDIR],data->end[JDIR]+joffset,
             data->beg[IDIR],data->end[IDIR]+ioffset,
    KOKKOS_LAMBDA (int k, int j, int i) {
      // Init the directions (should be in the kernel for proper optimisation by the compilers)
      EXPAND( const int Xn = DIR+MX1;                    ,
              const int Xt = (DIR == IDIR ? MX2 : MX1);  ,
              const int Xb = (DIR == KDIR ? MX2 : MX3);  )
      // Primitive variables
      real vL[Phys::nvar];
      real vR[Phys::nvar];
      real dv[Phys::nvar];

      // Conservative variables
      real uL[Phys::nvar];
      real uR[Phys::nvar];

      // Flux (left and right)
      real fluxL[Phys::nvar];
      real fluxR[Phys::nvar];

      // Roe
      real Rc[Phys::nvar][Phys::nvar];
      real um[Phys::nvar];

      // 1-- Store the primitive variables on the left, right, and averaged states
      extrapol.ExtrapolatePrimVar(i, j, k, vL, vR);
#pragma unroll
      for(int nv = 0 ; nv < Phys::nvar; nv++) {
        dv[nv] = vR[nv] - vL[nv];
      }

      // --- Compute the square of the sound speed
      real a, a2, a2L, a2R;
#if HAVE_ENERGY
      a2L = std::sqrt(eos.GetGamma(vL[PRS],vL[RHO])*(vL[PRS]/vL[RHO]));
      a2R = std::sqrt(eos.GetGamma(vR[PRS],vR[RHO])*(vR[PRS]/vR[RHO]));
      real h, vel2;
#else
      a2L = HALF_F*(eos.GetWaveSpeed(k,j,i)
                    +eos.GetWaveSpeed(k-koffset,j-joffset,i-ioffset));
      a2R = a2L;
#endif
      // Take the square
      a2L = a2L*a2L;
      a2R = a2R*a2R;

      // 2-- Compute the conservative variables
      K_PrimToCons<Phys>(uL, vL, &eos);
      K_PrimToCons<Phys>(uR, vR, &eos);

      // 3-- Compute the left and right fluxes
      K_Flux<Phys,DIR>(fluxL, vL, uL, a2L);
      K_Flux<Phys,DIR>(fluxR, vR, uR, a2R);

      // Compute gamma of this interface
      // todo(glesur): check that it's not the internal energy that should be used there instead
      #if HAVE_ENERGY
      real gamma = eos.GetGamma(0.5*(vL[PRS]+vR[PRS]), 0.5*(vL[RHO]+vR[RHO]));
      real gamma_m1 = gamma-1;
      #endif

      //  ----  Define Wave Jumps  ----
#if ROE_AVERAGE == YES
      real s, c;
      s       = std::sqrt(vR[RHO]/vL[RHO]);
      um[RHO] = vL[RHO]*s;
      s       = ONE_F/(ONE_F + s);
      c       = ONE_F - s;

      EXPAND(um[VX1] = s*vL[VX1] + c*vR[VX1];  ,
      um[VX2] = s*vL[VX2] + c*vR[VX2];  ,
      um[VX3] = s*vL[VX3] + c*vR[VX3];)

  #if HAVE_ENERGY
      real gmm1_inv = ONE_F / gamma_m1;

      vel2 = EXPAND(um[VX1]*um[VX1], + um[VX2]*um[VX2], + um[VX3]*um[VX3]);

      real hl, hr;
      hl  = HALF_F*(EXPAND(vL[VX1]*vL[VX1], + vL[VX2]*vL[VX2], + vL[VX3]*vL[VX3]));
      hl += a2L*gmm1_inv;

      hr = HALF_F*(EXPAND(vR[VX1]*vR[VX1], + vR[VX2]*vR[VX2], + vR[VX3]*vR[VX3]));
      hr += a2R*gmm1_inv;

      h = s*hl + c*hr;

      /* -------------------------------------------------
      the following should be  equivalent to

      scrh = EXPAND(   dv[VX1]*dv[VX1],
      + dv[VX2]*dv[VX2],
      + dv[VX3]*dv[VX3]);

      a2 = s*a2L + c*a2R + 0.5*gamma_m1*s*c*scrh;

      and therefore always positive.
      just work out the coefficiendnts...
      -------------------------------------------------- */

      a2 = gamma_m1*(h - HALF_F*vel2);
      a  = std::sqrt(a2);
  #else
      a2 = HALF_F*(a2L + a2R);
      a  = std::sqrt(a2);
  #endif // HAVE_ENERGY
#else
#pragma unroll
      for(int nv = 0 ; nv < Phys::nvar; nv++) {
        um[nv] = HALF_F*(vR[nv]+vL[nv]);
      }
  #if HAVE_ENERGY
      a2   = gamma*um[PRS]/um[RHO];
      a    = std::sqrt(a2);

      vel2 = EXPAND(um[VX1]*um[VX1], + um[VX2]*um[VX2], + um[VX3]*um[VX3]);
      h    = HALF_F*vel2 + a2/gamma_m1;
  #else
      a2 = HALF_F*(a2L + a2R);
      a  = std::sqrt(a2);
  #endif // HAVE_ENERGY
#endif // ROE_AVERAGE == YES/NO

// **********************************************************************************
      /* ----------------------------------------------------------------
      define non-zero components of conservative eigenvectors Rc,
      eigenvalues (lambda) and wave strenght eta = L.du
      ----------------------------------------------------------------  */

      real lambda[NMODES], alambda[NMODES];
      real eta[NMODES];

#pragma unroll
      for(int nv1 = 0 ; nv1 < Phys::nvar; nv1++) {
#pragma unroll
        for(int nv2 = 0 ; nv2 < Phys::nvar; nv2++) {
          Rc[nv1][nv2] = 0;
        }
      }

      //  ---- (u - c_s)  ----

      // nn         = 0;
      lambda[I0] = um[Xn] - a;
#if HAVE_ENERGY
      eta[I0] = HALF_F/a2*(dv[PRS] - dv[Xn]*um[RHO]*a);
#else
      eta[I0] = HALF_F*(dv[RHO] - um[RHO]*dv[Xn]/a);
#endif

      Rc[RHO][I0]        = ONE_F;

      EXPAND(Rc[Xn][I0] = um[Xn] - a;   ,
      Rc[Xt][I0] = um[Xt];       ,
      Rc[Xb][I0] = um[Xb];  )
#if HAVE_ENERGY
      Rc[ENG][I0] = h - um[Xn]*a;
#endif

      /*  ---- (u + c_s)  ----  */

      // nn         = 1;
      lambda[I1] = um[Xn] + a;
#if HAVE_ENERGY
      eta[I1]    = HALF_F/a2*(dv[PRS] + dv[Xn]*um[RHO]*a);
#else
      eta[I1] = HALF_F*(dv[RHO] + um[RHO]*dv[Xn]/a);
#endif

      Rc[RHO][I1]        = ONE_F;
      EXPAND(Rc[Xn][I1] = um[Xn] + a;   ,
      Rc[Xt][I1] = um[Xt];       ,
      Rc[Xb][I1] = um[Xb];)
#if HAVE_ENERGY
      Rc[ENG][I1] = h + um[Xn]*a;
#endif

#if HAVE_ENERGY
      /*  ----  (u)  ----  */

      // nn         = 2;
      lambda[IE] = um[Xn];
      eta[IE]    = dv[RHO] - dv[PRS]/a2;
      Rc[RHO][IE]        = ONE_F;
      EXPAND(Rc[MX1][IE] = um[VX1];   ,
      Rc[MX2][IE] = um[VX2];   ,
      Rc[MX3][IE] = um[VX3];)
      Rc[ENG][IE]        = HALF_F*vel2;
#endif

#if COMPONENTS > 1

      /*  ----  (u)  ----  */

      // nn++;
      lambda[I2] = um[Xn];
      eta[I2]    = um[RHO]*dv[Xt];
      Rc[Xt][I2] = ONE_F;
  #if HAVE_ENERGY
      Rc[ENG][I2] = um[Xt];
  #endif
#endif

#if COMPONENTS > 2

      /*  ----  (u)  ----  */

      // nn++;
      lambda[I3] = um[Xn];
      eta[I3]    = um[RHO]*dv[Xb];
      Rc[Xb][I3] = ONE_F;
  #if HAVE_ENERGY
      Rc[ENG][I3] = um[Xb];
  #endif
#endif

      /*  ----  get max eigenvalue  ----  */

      real cmax = FABS(um[Xn]) + a;
      //g_maxMach = FMAX(FABS(um[Xn]/a), g_maxMach);

      /* ---------------------------------------------
      use the HLL flux function if the interface
      lies within a strong shock.
      The effect of this switch is visible
      in the Mach reflection test.
      --------------------------------------------- */

      real scrh;
#if HAVE_ENERGY
      scrh  = FABS(vL[PRS] - vR[PRS]);
      scrh /= FMIN(vL[PRS],vR[PRS]);
#else
      scrh  = FABS(vL[RHO] - vR[RHO]);
      scrh /= FMIN(vL[RHO],vR[RHO]);
      scrh *= a*a;
#endif

/*#if CHECK_ROE_MATRIX == YES
      for(int nv = 0 ; nv < Phys::nvar; nv++) {
          um[nv] = ZERO_F;
          for(int nv1 = 0 ; nv1 < Phys::nvar; nv1++) {
              for(int nv2 = 0 ; nv2 < Phys::nvar; nv2++) {
                  um[nv] += Rc[nv][k]*(k==j)*lambda[k]*eta[j];
              }
          }
      }
      for(int nv = 0 ; nv < Phys::nvar; nv++) {
          scrh = fluxR[nv] - fluxL[nv] - um[nv];
          if (nv == Xn) scrh += pR - pL;
          if (FABS(scrh) > 1.e-6){
              print ("! Matrix condition not satisfied %d, %12.6e\n", nv, scrh);
              exit(1);
          }
      }
#endif*/

      if (scrh > HALF_F && (vR[Xn] < vL[Xn])) {   /* -- tunable parameter -- */
#if DIMENSIONS > 1
        real scrh1;
        real bmin, bmax;
        bmin = FMIN(ZERO_F, lambda[0]);
        bmax = FMAX(ZERO_F, lambda[1]);
        scrh1 = ONE_F/(bmax - bmin);
#pragma unroll
        for(int nv = 0 ; nv < Phys::nvar; nv++) {
          Flux(nv,k,j,i)  = bmin*bmax*(uR[nv] - uL[nv])
                  +   bmax*fluxL[nv] - bmin*fluxR[nv];
          Flux(nv,k,j,i) *= scrh1;
        }
#endif
      } else {
        /* -----------------------------------------------------------
                            compute Roe flux
        ----------------------------------------------------------- */

#pragma unroll
        for(int nv = 0 ; nv < Phys::nvar; nv++) {
          alambda[nv]  = fabs(lambda[nv]);
        }

        /*  ----  entropy fix  ----  */
        real delta = 1.e-7;
        if (alambda[0] <= delta) {
          alambda[0] = HALF_F*lambda[0]*lambda[0]/delta + HALF_F*delta;
        }
        if (alambda[1] <= delta) {
          alambda[1] = HALF_F*lambda[1]*lambda[1]/delta + HALF_F*delta;
        }

#pragma unroll
        for(int nv = 0 ; nv < Phys::nvar; nv++) {
          Flux(nv,k,j,i) = fluxL[nv] + fluxR[nv];
#pragma unroll
          for(int nv2 = 0 ; nv2 < Phys::nvar; nv2++) {
            Flux(nv,k,j,i) -= alambda[nv2]*eta[nv2]*Rc[nv][nv2];
          }
          Flux(nv,k,j,i) *= HALF_F;
        }
      }

      //6-- Compute maximum wave speed for this sweep
      cMax(k,j,i) = cmax;
    }
  );

  idfx::popRegion();
}

#endif // FLUID_RIEMANNSOLVER_HDSOLVERS_ROEHD_HPP_
