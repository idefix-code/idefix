// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_RIEMANNSOLVER_MHDSOLVERS_ROEMHD_HPP_
#define FLUID_RIEMANNSOLVER_MHDSOLVERS_ROEMHD_HPP_

#include "../idefix.hpp"
#include "extrapolateToFaces.hpp"
#include "flux.hpp"
#include "convertConsToPrim.hpp"
#include "storeFlux.hpp"
#include "constrainedTransport.hpp"

#define ROE_AVERAGE 0
#undef NMODES

#define KFASTM 0
#define KFASTP 1
#define KDIVB  2

#if HAVE_ENERGY

#define KENTRP 3
#define KSLOWM 4
#define KSLOWP 5
#define KALFVM 6
#define KALFVP 7

#define NMODES 8

#else

#define KSLOWM 3
#define KSLOWP 4
#define KALFVM 5
#define KALFVP 6

#define NMODES 7

#endif


/*! Return the sign of x. */
#define DSIGN(x) ( (x) >= 0.0 ? (1.0) : (-1.0))

// Compute Riemann fluxes from states using ROE solver
template <typename Phys>
template<const int DIR>
void RiemannSolver<Phys>::RoeMHD(IdefixArray4D<real> &Flux) {
  idfx::pushRegion("RiemannSolver::ROE_MHD");

  using EMF = ConstrainedTransport<Phys>;

  constexpr int ioffset = (DIR==IDIR) ? 1 : 0;
  constexpr int joffset = (DIR==JDIR) ? 1 : 0;
  constexpr int koffset = (DIR==KDIR) ? 1 : 0;

  // extension in perp to the direction of integration, as required by CT.
  constexpr int iextend = (DIR==IDIR) ? 0 : 1;
  #if DIMENSIONS > 1
    constexpr int jextend = (DIR==JDIR) ? 0 : 1;
  #else
    constexpr int jextend = 0;
  #endif
  #if DIMENSIONS > 2
    constexpr int kextend = (DIR==KDIR) ? 0 : 1;
  #else
    constexpr int kextend = 0;
  #endif

  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Vs = this->Vs;
  IdefixArray3D<real> cMax = this->cMax;

  // Required for high order interpolations
  IdefixArray1D<real> dx = this->data->dx[DIR];

  // References to required emf components
  IdefixArray3D<real> Eb;
  IdefixArray3D<real> Et;

  const typename EMF::AveragingType emfAverage = hydro->emf->averaging;


  // Required by UCT_Contact
  IdefixArray3D<real> SV;

  // Required by UCT_HLLX
  IdefixArray3D<real> aL;
  IdefixArray3D<real> aR;
  IdefixArray3D<real> dL;
  IdefixArray3D<real> dR;

  EquationOfState eos = *(hydro->eos.get());

  // TODO(baghdads) what is this delta?
  real delta    = 1.e-6;

  // Define normal, tangent and bi-tanget indices
  // st and sb will be useful only when Hall is included
  real st = ONE_F, sb = ONE_F;

  ExtrapolateToFaces<Phys,DIR> extrapol = *this->GetExtrapolator<DIR>();

  switch(DIR) {
    case(IDIR):
      D_EXPAND(
                st = -ONE_F;  ,
                  ,

                sb = +ONE_F;  )

      Et = hydro->emf->ezi;
      Eb = hydro->emf->eyi;

      SV = hydro->emf->svx;

      aL = hydro->emf->axL;
      aR = hydro->emf->axR;

      dL = hydro->emf->dxL;
      dR = hydro->emf->dxR;

      break;
#if DIMENSIONS >= 2
    case(JDIR):
      D_EXPAND(
                st = +ONE_F;  ,
                              ,

                sb = -ONE_F;  )

      Et = hydro->emf->ezj;
      Eb = hydro->emf->exj;

      SV = hydro->emf->svy;

      aL = hydro->emf->ayL;
      aR = hydro->emf->ayR;

      dL = hydro->emf->dyL;
      dR = hydro->emf->dyR;

      break;
#endif
#if DIMENSIONS == 3
    case(KDIR):

      D_EXPAND(

                st = -ONE_F;  ,
                  ,
                sb = +ONE_F;  )

      Et = hydro->emf->eyk;
      Eb = hydro->emf->exk;

      SV = hydro->emf->svz;

      aL = hydro->emf->azL;
      aR = hydro->emf->azR;

      dL = hydro->emf->dzL;
      dR = hydro->emf->dzR;
      break;
#endif
    default:
      IDEFIX_ERROR("Wrong direction");
  }

  idefix_for("CalcRiemannFlux",
             data->beg[KDIR]-kextend,data->end[KDIR]+koffset+kextend,
             data->beg[JDIR]-jextend,data->end[JDIR]+joffset+jextend,
             data->beg[IDIR]-iextend,data->end[IDIR]+ioffset+iextend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      // Init the directions (should be in the kernel for proper optimisation by the compilers)
      EXPAND( const int Xn = DIR+MX1;                    ,
              const int Xt = (DIR == IDIR ? MX2 : MX1);  ,
              const int Xb = (DIR == KDIR ? MX2 : MX3);  )

      EXPAND( const int BXn = DIR+BX1;                    ,
              const int BXt = (DIR == IDIR ? BX2 : BX1);  ,
              const int BXb = (DIR == KDIR ? BX2 : BX3);   )

      // Primitive variables
      real vL[Phys::nvar];
      real vR[Phys::nvar];
      real dV[Phys::nvar];

      // Conservative variables
      real uL[Phys::nvar];
      real uR[Phys::nvar];
      [[maybe_unused]] real dU[Phys::nvar];

      // Flux (left and right)
      real fluxL[Phys::nvar];
      real fluxR[Phys::nvar];

      // Roe
      real Rc[Phys::nvar][Phys::nvar];


      // 1-- Store the primitive variables on the left, right, and averaged states
      extrapol.ExtrapolatePrimVar(i, j, k, vL, vR);
      vL[BXn] = Vs(DIR,k,j,i);
      vR[BXn] = vL[BXn];

#pragma unroll
      for(int nv = 0 ; nv < Phys::nvar; nv++) {
        dV[nv] = vR[nv] - vL[nv];
      }

      // 2-- Compute the conservative variables
      K_PrimToCons<Phys>(uL, vL, &eos);
      K_PrimToCons<Phys>(uR, vR, &eos);

      // --- Compute the square of the sound speed
      real a, a2, a2L, a2R;
      #if HAVE_ENERGY
        // These are actually not used, but are initialised to avoid warnings
        a2L = ONE_F;
        a2R = ONE_F;
        real gamma = eos.GetGamma(0.5*(vL[RHO]+vR[RHO]),0.5*(vL[PRS]+vR[PRS]));
      #else
        a2L = HALF_F*(eos.GetWaveSpeed(k,j,i)
                    +eos.GetWaveSpeed(k-koffset,j-joffset,i-ioffset));
        a2L = a2L*a2L;
        a2R = a2L;
      #endif

      // 3-- Compute the left and right fluxes
#pragma unroll
      for(int nv = 0 ; nv < Phys::nvar; nv++) {
        fluxL[nv] = uL[nv];
        fluxR[nv] = uR[nv];
        dU[nv] = uR[nv] - uL[nv];
      }
      K_Flux<Phys,DIR>(fluxL, vL, fluxL, a2L);
      K_Flux<Phys,DIR>(fluxR, vR, fluxR, a2R);

      // 5. Set eigenvectors components Rc = 0 initially
#pragma unroll
      for(int nv1 = 0 ; nv1 < Phys::nvar; nv1++) {
#pragma unroll
        for(int nv2 = 0 ; nv2 < Phys::nvar; nv2++) {
          Rc[nv1][nv2] = 0;
        }
      }

      real sqr_rho_L, sqr_rho_R, sl, sr, rho, sqrt_rho;

      // 6c. Compute Roe averages
      sqr_rho_L = std::sqrt(vL[RHO]);
      sqr_rho_R = std::sqrt(vR[RHO]);

      sl = sqr_rho_L/(sqr_rho_L + sqr_rho_R);
      sr = sqr_rho_R/(sqr_rho_L + sqr_rho_R);

      // sl = sr = 0.5;

      rho = sr*vL[RHO] + sl*vR[RHO];

      sqrt_rho = std::sqrt(rho);

      [[maybe_unused]] real u, v, w, Bx, By, Bz, sBx, bx, by, bz, bt2, b2, Btmag;

      EXPAND ( u = sl*vL[Xn] + sr*vR[Xn];  ,
               v = sl*vL[Xt] + sr*vR[Xt];  ,
               w = sl*vL[Xb] + sr*vR[Xb];  )

      EXPAND ( Bx = sr*vL[BXn] + sl*vR[BXn];  ,
               By = sr*vL[BXt] + sl*vR[BXt];  ,
               Bz = sr*vL[BXb] + sl*vR[BXb];  )

      sBx = (Bx >= 0.0 ? 1.0 : -1.0);

      EXPAND( bx = Bx/sqrt_rho;  ,
              by = By/sqrt_rho;  ,
              bz = Bz/sqrt_rho;  )

      bt2   = EXPAND(0.0  , + by*by, + bz*bz);
      b2    = bx*bx + bt2;
      Btmag = std::sqrt(bt2*rho);

      real X  = EXPAND(dV[BXn]*dV[BXn], + dV[BXt]*dV[BXt], + dV[BXb]*dV[BXb]);
      X /= (sqr_rho_L + sqr_rho_R)*(sqr_rho_L + sqr_rho_R)*2.0;


      [[maybe_unused]] real Bmag2L, Bmag2R, pL, pR;
      Bmag2L = EXPAND(vL[BX1]*vL[BX1] , + vL[BX2]*vL[BX2], + vL[BX3]*vL[BX3]);
      Bmag2R = EXPAND(vR[BX1]*vR[BX1] , + vR[BX2]*vR[BX2], + vR[BX3]*vR[BX3]);
#if HAVE_ENERGY
      pL  = vL[PRS] + HALF_F*Bmag2L;
      pR  = vR[PRS] + HALF_F*Bmag2R;
#else
      pL  = a2L*vL[RHO] + HALF_F*Bmag2L;
      pR  = a2R*vR[RHO] + HALF_F*Bmag2R;
#endif

      // 6d. Compute enthalpy and sound speed.
#if HAVE_ENERGY
      real vel2, HL, HR, H, Hgas;
      real vdm, BdB;

      vdm = EXPAND(u*dU[Xn],  + v*dU[Xt],  + w*dU[Xb]);
      BdB = EXPAND(Bx*dU[BXn], + By*dU[BXt], + Bz*dU[BXb]);

      vel2    = EXPAND(u*u, + v*v, + w*w);
      dV[PRS] = (gamma-1.0)*((0.5*vel2 - X)*dV[RHO] - vdm + dU[ENG] - BdB);

      HL   = (uL[ENG] + pL)/vL[RHO];
      HR   = (uR[ENG] + pR)/vR[RHO];
      H    = sl*HL + sr*HR;   // total enthalpy

      Hgas = H - b2;         // gas enthalpy

      a2 = (2.0 - gamma)*X + (gamma-1.0)*(Hgas - 0.5*vel2);
      if (a2 < 0.0) {
          //IDEFIX_ERROR("! Roe_Solver(): a2 < 0.0 !! \n");
      }
#else
      // in most cases a2L = a2R for isothermal MHD
      a2 = 0.5*(a2L + a2R) + X;
#endif

      /* ------------------------------------------------------------
      6e. Compute fast and slow magnetosonic speeds.

      The following expression appearing in the definitions
      of the fast magnetosonic speed

      (a^2 - b^2)^2 + 4*a^2*bt^2 = (a^2 + b^2)^2 - 4*a^2*bx^2

      is always positive and avoids round-off errors.

      Note that we always use the total field to compute the
      characteristic speeds.
      ------------------------------------------------------------ */

      [[maybe_unused]] real scrh, ca, cf, cs, ca2, cf2, cs2, alpha_f, alpha_s, beta_y, beta_z;
      scrh = a2 - b2;
      ca2  = bx*bx;
      scrh = scrh*scrh + 4.0*bt2*a2;
      scrh = std::sqrt(scrh);

      cf2 = 0.5*(a2 + b2 + scrh);
      cs2 = a2*ca2/cf2;   // -- same as 0.5*(a2 + b2 - scrh)

      cf = std::sqrt(cf2);
      cs = std::sqrt(cs2);
      ca = std::sqrt(ca2);
      a  = std::sqrt(a2);

      if (cf == cs) {
        alpha_f = 1.0;
        alpha_s = 0.0;
      } else if (a <= cs) {
        alpha_f = 0.0;
        alpha_s = 1.0;
      } else if (cf <= a) {
        alpha_f = 1.0;
        alpha_s = 0.0;
      } else {
        scrh    = 1.0/(cf2 - cs2);
        alpha_f = (a2  - cs2)*scrh;
        alpha_s = (cf2 -  a2)*scrh;
        alpha_f = FMAX(0.0, alpha_f);
        alpha_s = FMAX(0.0, alpha_s);
        alpha_f = std::sqrt(alpha_f);
        alpha_s = std::sqrt(alpha_s);
      }

      if (Btmag > 1.e-9) {
        SELECT(                      ,
                beta_y = DSIGN(By);  ,
                beta_y = By/Btmag;
                beta_z = Bz/Btmag;   )
      } else {
        SELECT(                         ,
                beta_y = 1.0;           ,
                beta_z = beta_y = 1.0;  )
      }

      /* -------------------------------------------------------------------
      6f. Compute non-zero entries of conservative eigenvectors (Rc),
          wave strength L*dU (=eta) for all 8 (or 7) waves using the
          expressions given by Eq. [4.18]--[4.21].
          Fast and slow eigenvectors are multiplied by a^2 while
          jumps are divided by a^2.

          Notes:
          - the expression on the paper has a typo in the very last term
          of the energy component: it should be + and not - !
          - with background field splitting: additional terms must be
          added to the energy component for fast, slow and Alfven waves.
          To obtain energy element, conservative eigenvector (with
          total field) must be multiplied by | 0 0 0 0 -B0y -B0z 1 |.
          Also, H - b2 does not give gas enthalpy. A term b0*btot must
          be added and eta (wave strength) should contain total field
          and deviation's delta.
      ------------------------------------------------------------------- */

      // Fast wave:  u - c_f
      real lambda[NMODES], alambda[NMODES], eta[NMODES];
      [[maybe_unused]] real beta_dv, beta_dB, beta_v;

      int kk = KFASTM;
      lambda[kk] = u - cf;

      scrh    = alpha_s*cs*sBx;
      beta_dv = EXPAND(0.0, + beta_y*dV[Xt], + beta_z*dV[Xb]);
      beta_dB = EXPAND(0.0, + beta_y*dV[BXt], + beta_z*dV[BXb]);

      Rc[RHO][kk] = alpha_f;
      EXPAND( Rc[Xn][kk] = alpha_f*lambda[kk];       ,
              Rc[Xt][kk] = alpha_f*v + scrh*beta_y;  ,
              Rc[Xb][kk] = alpha_f*w + scrh*beta_z;  )

      EXPAND(                                           ,
              Rc[BXt][kk] = alpha_s*a*beta_y/sqrt_rho;  ,
              Rc[BXb][kk] = alpha_s*a*beta_z/sqrt_rho;  )

#if HAVE_ENERGY
      beta_v  = EXPAND(0.0, + beta_y*v,       + beta_z*w);
      Rc[ENG][kk] =   alpha_f*(Hgas - u*cf) + scrh*beta_v
                  + alpha_s*a*Btmag/sqrt_rho;

      eta[kk] =   alpha_f*(X*dV[RHO] + dV[PRS]);
#else
      // eta[kk] =   alpha_f*(0.0*X + a2)*dV[RHO] + rho*scrh*beta_dv
      //        - rho*alpha_f*cf*dV[Xn] + sqrt_rho*alpha_s*a*beta_dB;
      eta[kk] =   alpha_f*a2*dV[RHO];
#endif
      eta[kk] += rho*scrh*beta_dv - rho*alpha_f*cf*dV[Xn] + sqrt_rho*alpha_s*a*beta_dB;
      eta[kk] *= 0.5/a2;

      // Fast wave:  u + c_f

      kk = KFASTP;
      lambda[kk] = u + cf;

      Rc[RHO][kk] = alpha_f;
      EXPAND( Rc[Xn][kk] = alpha_f*lambda[kk];       ,
              Rc[Xt][kk] = alpha_f*v - scrh*beta_y;  ,
              Rc[Xb][kk] = alpha_f*w - scrh*beta_z;  )
      EXPAND(                                 ,
              Rc[BXt][kk] = Rc[BXt][KFASTM];  ,
              Rc[BXb][kk] = Rc[BXb][KFASTM];  )

#if HAVE_ENERGY
      Rc[ENG][kk] =   alpha_f*(Hgas + u*cf) - scrh*beta_v
                  + alpha_s*a*Btmag/sqrt_rho;

      eta[kk] =   alpha_f*(X*dV[RHO] + dV[PRS]) - rho*scrh*beta_dv
              + rho*alpha_f*cf*dV[Xn]        + sqrt_rho*alpha_s*a*beta_dB;
#else
      eta[kk] =   alpha_f*(0.*X + a2)*dV[RHO] - rho*scrh*beta_dv
              + rho*alpha_f*cf*dV[Xn]      + sqrt_rho*alpha_s*a*beta_dB;
#endif

      eta[kk] *= 0.5/a2;

      // Entropy wave:  u

#if HAVE_ENERGY
      kk = KENTRP;
      lambda[kk] = u;

      Rc[RHO][kk] = 1.0;
      EXPAND( Rc[Xn][kk] = u;  ,
              Rc[Xt][kk] = v;  ,
              Rc[Xb][kk] = w;  )
      Rc[ENG][kk] = 0.5*vel2 + (gamma - 2.0)/(gamma-1.0)*X;

      eta[kk] = ((a2 - X)*dV[RHO] - dV[PRS])/a2;
#endif

      /* -----------------------------------------------------------------
      div.B wave (u): this wave exists when:

      1) 8 wave formulation
      2) CT, since we always have 8 components, but it
          carries zero jump.

      With GLM, KDIVB is replaced by KPSI_GLMM, KPSI_GLMP and these
      two waves should not enter in the Riemann solver (eta = 0.0)
      since the 2x2 linear system formed by (B,psi) has already
      been solved.
      ----------------------------------------------------------------- */

      kk = KDIVB;
      lambda[kk] = u;
      eta[kk] = 0.0;

#if COMPONENTS > 1
      // Slow wave:  u - c_s

      scrh = alpha_f*cf*sBx;

      kk = KSLOWM;
      lambda[kk] = u - cs;

      Rc[RHO][kk] = alpha_s;
      EXPAND( Rc[Xn][kk] = alpha_s*lambda[kk];       ,
              Rc[Xt][kk] = alpha_s*v - scrh*beta_y;  ,
              Rc[Xb][kk] = alpha_s*w - scrh*beta_z;  )
      EXPAND(                                             ,
              Rc[BXt][kk] = - alpha_f*a*beta_y/sqrt_rho;  ,
              Rc[BXb][kk] = - alpha_f*a*beta_z/sqrt_rho;  )

  #if HAVE_ENERGY
      Rc[ENG][kk] =   alpha_s*(Hgas - u*cs) - scrh*beta_v
                  - alpha_f*a*Btmag/sqrt_rho;

      eta[kk] =   alpha_s*(X*dV[RHO] + dV[PRS]) - rho*scrh*beta_dv
              - rho*alpha_s*cs*dV[Xn]        - sqrt_rho*alpha_f*a*beta_dB;
  #else
      eta[kk] =   alpha_s*(0.*X + a2)*dV[RHO] - rho*scrh*beta_dv
              - rho*alpha_s*cs*dV[Xn]      - sqrt_rho*alpha_f*a*beta_dB;
  #endif

      eta[kk] *= 0.5/a2;

      // Slow wave:  u + c_s

      kk = KSLOWP;
      lambda[kk] = u + cs;

      Rc[RHO][kk] = alpha_s;
      EXPAND( Rc[Xn][kk] = alpha_s*lambda[kk];       ,
              Rc[Xt][kk] = alpha_s*v + scrh*beta_y;  ,
              Rc[Xb][kk] = alpha_s*w + scrh*beta_z;  )
      EXPAND(                                 ,
              Rc[BXt][kk] = Rc[BXt][KSLOWM];  ,
              Rc[BXb][kk] = Rc[BXb][KSLOWM];  )

  #if HAVE_ENERGY
      Rc[ENG][kk] =   alpha_s*(Hgas + u*cs) + scrh*beta_v
                  - alpha_f*a*Btmag/sqrt_rho;

      eta[kk] =   alpha_s*(X*dV[RHO] + dV[PRS]) + rho*scrh*beta_dv
              + rho*alpha_s*cs*dV[Xn]        - sqrt_rho*alpha_f*a*beta_dB;
  #else
      eta[kk] =   alpha_s*(0.*X + a2)*dV[RHO] + rho*scrh*beta_dv
              + rho*alpha_s*cs*dV[Xn]      - sqrt_rho*alpha_f*a*beta_dB;
  #endif

      eta[kk] *= 0.5/a2;

#endif // COMPONENTS > 1

#if COMPONENTS == 3

      // Alfven wave:  u - c_a

      kk = KALFVM;
      lambda[kk] = u - ca;

      Rc[Xt][kk] = - rho*beta_z;
      Rc[Xb][kk] = + rho*beta_y;
      Rc[BXt][kk] = - sBx*sqrt_rho*beta_z;
      Rc[BXb][kk] =   sBx*sqrt_rho*beta_y;
  #if HAVE_ENERGY
      Rc[ENG][kk] = - rho*(v*beta_z - w*beta_y);
  #endif

      eta[kk] = + beta_y*dV[Xb]               - beta_z*dV[Xt]
              + sBx/sqrt_rho*(beta_y*dV[BXb] - beta_z*dV[BXt]);

      eta[kk] *= 0.5;

      // Alfven wave:  u + c_a

      kk = KALFVP;
      lambda[kk] = u + ca;

      Rc[Xt][kk] = - Rc[Xt][KALFVM];
      Rc[Xb][kk] = - Rc[Xb][KALFVM];
      Rc[BXt][kk] =   Rc[BXt][KALFVM];
      Rc[BXb][kk] =   Rc[BXb][KALFVM];
  #if HAVE_ENERGY
      Rc[ENG][kk] = - Rc[ENG][KALFVM];
  #endif

      eta[kk] = - beta_y*dV[Xb]               + beta_z*dV[Xt]
              + sBx/sqrt_rho*(beta_y*dV[BXb] - beta_z*dV[BXt]);

      eta[kk] *= 0.5;
#endif // COMPONENTS == 3

      // 6g. Compute maximum signal velocity

      real cmax = std::fabs(u) + cf;

      // 6h. Save max and min Riemann fan speeds for EMF computation.
      sl = lambda[KFASTM];
      sr = lambda[KFASTP];

#pragma unroll
      for(int nv = 0 ; nv < Phys::nvar; nv++) {
          alambda[nv] = fabs(lambda[nv]);
      }

      // 6i. Entropy Fix

      if (alambda[KFASTM] < 0.5*delta) {
        alambda[KFASTM] = lambda[KFASTM]*lambda[KFASTM]/delta + 0.25*delta;
      }
      if (alambda[KFASTP] < 0.5*delta) {
        alambda[KFASTP] = lambda[KFASTP]*lambda[KFASTP]/delta + 0.25*delta;
      }
#if COMPONENTS > 1
      if (alambda[KSLOWM] < 0.5*delta) {
        alambda[KSLOWM] = lambda[KSLOWM]*lambda[KSLOWM]/delta + 0.25*delta;
      }
      if (alambda[KSLOWP] < 0.5*delta) {
        alambda[KSLOWP] = lambda[KSLOWP]*lambda[KSLOWP]/delta + 0.25*delta;
      }
#endif

      // 6j. Compute Roe numerical flux
#pragma unroll
      for(int nv1 = 0 ; nv1 < Phys::nvar; nv1++) {
        scrh = 0.0;
#pragma unroll
        for(int nv2 = 0 ; nv2 < Phys::nvar; nv2++) {
          scrh += alambda[nv2]*eta[nv2]*Rc[nv1][nv2];
        }
        Flux(nv1,k,j,i) = 0.5*(fluxL[nv1] + fluxR[nv1] - scrh);
      }

      // save maximum wave speed for this sweep
      cMax(k,j,i) = cmax;

      // 7-- Store the flux in the emf components
      if (emfAverage==EMF::arithmetic
                || emfAverage==EMF::uct0) {
        K_StoreEMF<DIR>(i,j,k,st,sb,Flux,Et,Eb);
      } else if (emfAverage==EMF::uct_contact) {
        K_StoreContact<DIR>(i,j,k,st,sb,Flux,Et,Eb,SV);
      } else if (emfAverage==EMF::uct_hll) {
        K_StoreHLL<DIR>(i,j,k,st,sb,sl,sr,vL,vR,Et,Eb,aL,aR,dL,dR);
      } else if (emfAverage==EMF::uct_hlld) {
        K_StoreHLLD<DIR>(i,j,k,st,sb,a2L,sl,sr,
                         vL,vR,uL,uR,Et,Eb,aL,aR,dL,dR);
      }
  });

  idfx::popRegion();
}
#endif // FLUID_RIEMANNSOLVER_MHDSOLVERS_ROEMHD_HPP_
