// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_HDSOLVERS_ROEHD_HPP_
#define HYDRO_HDSOLVERS_ROEHD_HPP_

#include "../idefix.hpp"
#include "hydro.hpp"
#include "slopeLimiter.hpp"
#include "fluxHD.hpp"
#include "convertConsToPrimHD.hpp"

#define ROE_AVERAGE 0

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
template<const int DIR>
void Hydro::RoeHD() {
  idfx::pushRegion("Hydro::ROE_Solver");

  int ioffset,joffset,koffset;
  // Determine the offset along which we do the extrapolation
  switch(DIR) {
    case(IDIR):
      ioffset = 1;
      joffset=koffset=0;
      break;
    case(JDIR):
      joffset=1;
      ioffset=koffset=0;
      break;
    case(KDIR):
      koffset=1;
      ioffset=joffset=0;
      break;
    default:
      IDEFIX_ERROR("Wrong direction");
  }

  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Vs = this->Vs;
  IdefixArray4D<real> Flux = this->FluxRiemann;
  IdefixArray3D<real> cMax = this->cMax;
  IdefixArray3D<real> csIsoArr = this->isoSoundSpeedArray;

  [[maybe_unused]] real gamma = this->gamma;
  [[maybe_unused]] real gamma_m1 = this->gamma - ONE_F;
  [[maybe_unused]] real csIso = this->isoSoundSpeed;
  [[maybe_unused]] HydroModuleStatus haveIsoCs = this->haveIsoSoundSpeed;

  SlopeLimiter<DIR,NVAR> slopeLim(Vc,data->dx[DIR],shockFlattening);

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
      real vL[NVAR];
      real vR[NVAR];
      real dv[NVAR];

      // Conservative variables
      real uL[NVAR];
      real uR[NVAR];

      // Flux (left and right)
      real fluxL[NVAR];
      real fluxR[NVAR];

      // Roe
      real Rc[NVAR][NVAR];
      real um[NVAR];

      // 1-- Store the primitive variables on the left, right, and averaged states
      slopeLim.ExtrapolatePrimVar(i, j, k, vL, vR);
#pragma unroll
      for(int nv = 0 ; nv < NVAR; nv++) {
        dv[nv] = vR[nv] - vL[nv];
      }

      // --- Compute the square of the sound speed
      real a, a2, a2L, a2R;
#if HAVE_ENERGY
      a2L = gamma * vL[PRS] / vL[RHO];
      a2R = gamma * vR[PRS] / vR[RHO];
      real h, vel2;
#else
      if(haveIsoCs == UserDefFunction) {
        a2L = HALF_F*(csIsoArr(k,j,i)+csIsoArr(k-koffset,j-joffset,i-ioffset));
      } else {
        a2L = csIso;
      }
      // Take the square
      a2L = a2L*a2L;
      a2R = a2L;
#endif

      // 2-- Compute the conservative variables
      K_PrimToCons(uL, vL, gamma_m1);
      K_PrimToCons(uR, vR, gamma_m1);

      // 3-- Compute the left and right fluxes
      K_Flux(fluxL, vL, uL, a2L, Xn);
      K_Flux(fluxR, vR, uR, a2R, Xn);

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
      for(int nv = 0 ; nv < NVAR; nv++) {
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
      for(int nv1 = 0 ; nv1 < NVAR; nv1++) {
#pragma unroll
        for(int nv2 = 0 ; nv2 < NVAR; nv2++) {
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
      for(int nv = 0 ; nv < NVAR; nv++) {
          um[nv] = ZERO_F;
          for(int nv1 = 0 ; nv1 < NVAR; nv1++) {
              for(int nv2 = 0 ; nv2 < NVAR; nv2++) {
                  um[nv] += Rc[nv][k]*(k==j)*lambda[k]*eta[j];
              }
          }
      }
      for(int nv = 0 ; nv < NVAR; nv++) {
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
        for(int nv = 0 ; nv < NVAR; nv++) {
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
        for(int nv = 0 ; nv < NVAR; nv++) {
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
        for(int nv = 0 ; nv < NVAR; nv++) {
          Flux(nv,k,j,i) = fluxL[nv] + fluxR[nv];
#pragma unroll
          for(int nv2 = 0 ; nv2 < NVAR; nv2++) {
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

#endif // HYDRO_HDSOLVERS_ROEHD_HPP_
