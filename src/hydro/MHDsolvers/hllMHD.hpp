// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_MHDSOLVERS_HLLMHD_HPP_
#define HYDRO_MHDSOLVERS_HLLMHD_HPP_

#include "../idefix.hpp"
#include "slopeLimiter.hpp"
#include "fluxMHD.hpp"
#include "convertConsToPrimMHD.hpp"
#include "storeFlux.hpp"
#include "electroMotiveForce.hpp"

// Compute Riemann fluxes from states using HLL solver
template<const int DIR>
void Hydro::HllMHD() {
  idfx::pushRegion("Hydro::HLL_MHD");

  int ioffset,joffset,koffset;
  int iextend, jextend,kextend;
  ioffset=joffset=koffset=0;
  // extension in perp to the direction of integration, as required by CT.
  iextend=jextend=kextend=0;

  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Vs = this->Vs;
  IdefixArray4D<real> Flux = this->FluxRiemann;
  IdefixArray3D<real> cMax = this->cMax;
  IdefixArray3D<real> csIsoArr = this->isoSoundSpeedArray;

  HydroModuleStatus haveHall = this->hallStatus.status;
  IdefixArray4D<real> J = this->J;
  IdefixArray3D<real> xHallArr = this->xHall;
  IdefixArray1D<real> dx = data->dx[DIR];
  IdefixArray1D<real> dx2 = data->dx[JDIR];
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> rt = data->rt;
  IdefixArray1D<real> dmu = data->dmu;

  // References to required emf components
  IdefixArray3D<real> Eb;
  IdefixArray3D<real> Et;

  const ElectroMotiveForce::AveragingType emfAverage = emf.averaging;

  // Required by UCT_Contact
  IdefixArray3D<int> SV;

  // Required by UCT_HLLX
  IdefixArray3D<real> aL;
  IdefixArray3D<real> aR;
  IdefixArray3D<real> dL;
  IdefixArray3D<real> dR;

  real gamma = this->gamma;
  [[maybe_unused]] real xHConstant = this->xH;
  [[maybe_unused]] real gamma_m1=gamma-ONE_F;
  [[maybe_unused]] real csIso = this->isoSoundSpeed;
  [[maybe_unused]] HydroModuleStatus haveIsoCs = this->haveIsoSoundSpeed;

  SlopeLimiter<DIR,NVAR> slopeLim(Vc,data->dx[DIR],shockFlattening);

  // Define normal, tangent and bi-tanget indices
  // st and sb will be useful only when Hall is included
  real st = ONE_F, sb = ONE_F;

  switch(DIR) {
    case(IDIR):
      ioffset = 1;
      D_EXPAND(
                st = -ONE_F;  ,
                jextend = 1;  ,
                kextend = 1;
                sb = +ONE_F;  )

      Et = this->emf.ezi;
      Eb = this->emf.eyi;

      SV = this->emf.svx;

      aL = this->emf.axL;
      aR = this->emf.axR;

      dL = this->emf.dxL;
      dR = this->emf.dxR;

      break;
#if DIMENSIONS >= 2
    case(JDIR):
      joffset=1;
      D_EXPAND(
                iextend = 1;
                st = +ONE_F;  ,
                              ,
                kextend = 1;
                sb = -ONE_F;  )

      Et = this->emf.ezj;
      Eb = this->emf.exj;

      SV = this->emf.svy;

      aL = this->emf.ayL;
      aR = this->emf.ayR;

      dL = this->emf.dyL;
      dR = this->emf.dyR;

      break;
#endif
#if DIMENSIONS == 3
    case(KDIR):
      koffset=1;
      D_EXPAND(
                iextend = 1;
                st = -ONE_F;  ,
                jextend = 1;  ,
                sb = +ONE_F;  )

      Et = this->emf.eyk;
      Eb = this->emf.exk;

      SV = this->emf.svz;

      aL = this->emf.azL;
      aR = this->emf.azR;

      dL = this->emf.dzL;
      dR = this->emf.dzR;
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
      real vL[NVAR];
      real vR[NVAR];

      // Conservative variables
      real uL[NVAR];
      real uR[NVAR];

      // Flux (left and right)
      real fluxL[NVAR];
      real fluxR[NVAR];

      // Signal speeds
      real cL, cR, cmax, c2Iso;

      c2Iso = ZERO_F;

      // 1-- Store the primitive variables on the left, right, and averaged states
      slopeLim.ExtrapolatePrimVar(i, j, k, vL, vR);
      vL[BXn] = Vs(DIR,k,j,i);
      vR[BXn] = vL[BXn];

      // 2-- Get the wave speed
      real gpr, b1, b2, b3, Btmag2, Bmag2;
      real xH;
#if HAVE_ENERGY
      gpr = gamma*vL[PRS];
#else
      if(haveIsoCs == UserDefFunction) {
        c2Iso = HALF_F*(csIsoArr(k,j,i)+csIsoArr(k-koffset,j-joffset,i-ioffset));
        c2Iso = c2Iso*c2Iso;
      } else {
        c2Iso = csIso*csIso;
      }

      gpr = c2Iso*vL[RHO];
#endif

      // -- get total field
      b1 = b2 = b3 = ZERO_F;
      EXPAND (b1 = vL[BXn];  ,
              b2 = vL[BXt];  ,
              b3 = vL[BXb];)

      Btmag2 = b2*b2 + b3*b3;
      Bmag2  = b1*b1 + Btmag2;

      cL = gpr - Bmag2;
      cL = gpr + Bmag2 + std::sqrt(cL*cL + FOUR_F*gpr*Btmag2);
      cL = std::sqrt(HALF_F*cL/vL[RHO]);

#if HAVE_ENERGY
      gpr = gamma*vR[PRS];
#else
      gpr = c2Iso*vR[RHO];
#endif

      // -- get total field
      b1 = b2 = b3 = ZERO_F;
      EXPAND (b1 = vR[BXn];  ,
              b2 = vR[BXt];  ,
              b3 = vR[BXb];)

      Btmag2 = b2*b2 + b3*b3;
      Bmag2  = b1*b1 + Btmag2;

      cR = gpr - Bmag2;
      cR = gpr + Bmag2 + std::sqrt(cR*cR + FOUR_F*gpr*Btmag2);
      cR = std::sqrt(HALF_F*cR/vR[RHO]);

      // 4.1
      real cminL = vL[Xn] - cL;
      real cmaxL = vL[Xn] + cL;

      real cminR = vR[Xn] - cR;
      real cmaxR = vR[Xn] + cR;

      real sl = FMIN(cminL, cminR);
      real sr = FMAX(cmaxL, cmaxR);

      // Signal speeds specific to B (different from the other ones when Hall is enabled)
      real SLb = sl;
      real SRb = sr;
      // if Hall is enabled, add whistler speed to the fan
      if(haveHall) {
        // Compute xHall
        if(haveHall==UserDefFunction) {
            if(DIR==IDIR) xH = AVERAGE_3D_X(xHallArr,k,j,i);
            if(DIR==JDIR) xH = AVERAGE_3D_Y(xHallArr,k,j,i);
            if(DIR==KDIR) xH = AVERAGE_3D_Z(xHallArr,k,j,i);
        } else {
          xH = xHConstant;
        }

        const int ig = ioffset*i + joffset*j + koffset*k;
        real dl = dx(ig);
        #if GEOMETRY == POLAR
            if(DIR==JDIR) dl = dl*x1(i);
        #elif GEOMETRY == SPHERICAL
            if(DIR==JDIR) dl = dl*rt(i);
            if(DIR==KDIR) dl = dl*rt(i)*dmu(j)/dx2(j);
        #endif

        real cw = FABS(xH) * std::sqrt(Bmag2) / dl;

        cminL = cminL - cw;
        cmaxL = cmaxL + cw;
        cminR = cminR - cw;
        cmaxR = cmaxR + cw;

        SLb = FMIN(cminL, cminR);
        SRb = FMAX(cmaxL,cmaxR);
      }

      cmax = FMAX(FABS(SLb), FABS(SRb));

      // 2-- Compute the conservative variables
      K_PrimToCons(uL, vL, gamma_m1);
      K_PrimToCons(uR, vR, gamma_m1);



#pragma unroll
      for(int nv = 0 ; nv < NVAR; nv++) {
        fluxL[nv] = uL[nv];
        fluxR[nv] = uR[nv];
      }



      // 3-- Compute the left and right fluxes
      K_Flux(fluxL, vL, fluxL, c2Iso, ARG_EXPAND(Xn, Xt, Xb), ARG_EXPAND(BXn, BXt, BXb));
      K_Flux(fluxR, vR, fluxR, c2Iso, ARG_EXPAND(Xn, Xt, Xb), ARG_EXPAND(BXn, BXt, BXb));

      // 4-- Compute the Hall flux
      if(haveHall) {
        [[maybe_unused]] int ip1, jp1, kp1;
        real Jx1, Jx2, Jx3;
        ip1=i+1;
        #if DIMENSIONS >=2
            jp1 = j+1;
        #else
            jp1=j;
        #endif
        #if DIMENSIONS == 3
            kp1 = k+1;
        #else
            kp1 = k;
        #endif

        if(DIR == IDIR) {
          Jx1 = AVERAGE_4D_XYZ(J, IDIR, kp1, jp1, i);
          Jx2 = AVERAGE_4D_Z(J, JDIR, kp1, j, i);
          Jx3 = AVERAGE_4D_Y(J, KDIR, k, jp1, i);
          #if COMPONENTS >= 2
          fluxL[BX2] += -xH* (  Jx1*uL[BX2] - Jx2*uL[BX1] );
          fluxR[BX2] += -xH* (  Jx1*uR[BX2] - Jx2*uR[BX1] );
          #endif
          #if COMPONENTS == 3
          fluxL[BX3] += -xH* (  Jx1*uL[BX3] - Jx3*uL[BX1] );
          fluxR[BX3] += -xH* (  Jx1*uR[BX3] - Jx3*uR[BX1] );
          #endif
        }
        if(DIR == JDIR) {
          Jx1 = AVERAGE_4D_Z(J, IDIR, kp1, j, i);
          Jx2 = AVERAGE_4D_XYZ(J, JDIR, kp1, j, ip1);
          Jx3 = AVERAGE_4D_X(J, KDIR, k, j, ip1);

          #if COMPONENTS >= 2
          fluxL[BX1] += -xH* (  Jx2*uL[BX1] - Jx1*uL[BX2] );
          fluxR[BX1] += -xH* (  Jx2*uR[BX1] - Jx1*uR[BX2] );
          #endif

          #if COMPONENTS == 3
          fluxL[BX3] += -xH* (  Jx2*uL[BX3] - Jx3*uL[BX2] );
          fluxR[BX3] += -xH* (  Jx2*uR[BX3] - Jx3*uR[BX2] );
          #endif
        }
        if(DIR == KDIR) {
          Jx1 = AVERAGE_4D_Y(J, IDIR, k, jp1, i);
          Jx2 = AVERAGE_4D_X(J, JDIR, k, j, ip1);
          Jx3 = AVERAGE_4D_XYZ(J, KDIR, k, jp1, ip1);

          fluxL[BX1] += -xH* (  Jx3*uL[BX1]  );
          fluxR[BX1] += -xH* (  Jx3*uR[BX1]  );

          #if COMPONENTS >= 2
          fluxL[BX2] += -xH* (  Jx3*uL[BX2]  );
          fluxR[BX2] += -xH* (  Jx3*uR[BX2]  );

          #if COMPONENTS==3
          fluxL[BX1] += -xH* (   - Jx1*uL[BX3] );
          fluxR[BX1] += -xH* (   - Jx1*uR[BX3] );

          fluxL[BX2] += -xH* (   - Jx2*uL[BX3] );
          fluxR[BX2] += -xH* (   - Jx2*uR[BX3] );
          #endif
          #endif
        }

        #if HAVE_ENERGY
          real JB = EXPAND(uL[BX1]*Jx1,  +uL[BX2]*Jx2, +uL[BX3]*Jx3 );
          real b2 = HALF_F*(EXPAND(uL[BX1]*uL[BX1], +uL[BX2]*uL[BX2], +uL[BX3]*uL[BX3]));
          if(DIR == IDIR) fluxL[ENG] += -xH* (Jx1*b2 - JB*uL[BX1]);
          #if COMPONENTS>=2
          if(DIR == JDIR) fluxL[ENG] += -xH* (Jx2*b2 - JB*uL[BX2]);
          #endif
          #if COMPONENTS >=3
          if(DIR == KDIR) fluxL[ENG] += -xH* (Jx3*b2 - JB*uL[BX3]);
          #endif

          JB = EXPAND(uR[BX1]*Jx1,  +uR[BX2]*Jx2, +uR[BX3]*Jx3 );
          b2 = HALF_F*(EXPAND(uR[BX1]*uR[BX1], +uR[BX2]*uR[BX2], +uR[BX3]*uR[BX3]));
          if(DIR == IDIR) fluxR[ENG] += -xH* (Jx1*b2 - JB*uR[BX1]);
          #if COMPONENTS>=2
          if(DIR == JDIR) fluxR[ENG] += -xH* (Jx2*b2 - JB*uR[BX2]);
          #endif
          #if COMPONENTS >=3
          if(DIR == KDIR) fluxR[ENG] += -xH* (Jx3*b2 - JB*uR[BX3]);
          #endif
      #endif
      }

      // 5-- Compute the flux from the left and right states
      if (sl > 0) {
        Flux(RHO,k,j,i) = fluxL[RHO];
        EXPAND( Flux(MX1,k,j,i) = fluxL[MX1];  ,
                Flux(MX2,k,j,i) = fluxL[MX2];  ,
                Flux(MX3,k,j,i) = fluxL[MX3];  )
      } else if (sr < 0) {
        Flux(RHO,k,j,i) = fluxR[RHO];
        EXPAND( Flux(MX1,k,j,i) = fluxR[MX1];  ,
                Flux(MX2,k,j,i) = fluxR[MX2];  ,
                Flux(MX3,k,j,i) = fluxR[MX3];  )
      } else {
        Flux(RHO,k,j,i) = (sl*sr*uR[RHO] - sl*sr*uL[RHO] + sr*fluxL[RHO] - sl*fluxR[RHO])
                          / (sr - sl);
        EXPAND( Flux(MX1,k,j,i) = (sl*sr*uR[MX1] - sl*sr*uL[MX1] + sr*fluxL[MX1] - sl*fluxR[MX1])
                                  / (sr - sl);  ,
                Flux(MX2,k,j,i) = (sl*sr*uR[MX2] - sl*sr*uL[MX2] + sr*fluxL[MX2] - sl*fluxR[MX2])
                                  / (sr - sl);  ,
                Flux(MX3,k,j,i) = (sl*sr*uR[MX3] - sl*sr*uL[MX3] + sr*fluxL[MX3] - sl*fluxR[MX3])
                                  / (sr - sl);  )
      }

      if (SLb > 0) {
#pragma unroll
        for (int nv = BX1 ; nv < NFLX; nv++) {
          Flux(nv,k,j,i) = fluxL[nv];
        }
      } else if (SRb < 0) {
#pragma unroll
        for (int nv = BX1 ; nv < NFLX; nv++) {
          Flux(nv,k,j,i) = fluxR[nv];
        }
      } else {
#pragma unroll
        for(int nv = BX1 ; nv < NFLX; nv++) {
          Flux(nv,k,j,i) = SLb*SRb*uR[nv] - SLb*SRb*uL[nv] + SRb*fluxL[nv] - SLb*fluxR[nv];
          Flux(nv,k,j,i) *= (1.0 / (SRb - SLb));
        }
      }


      //6-- Compute maximum wave speed for this sweep
      cMax(k,j,i) = cmax;

      // 7-- Store the flux in the emf components
      if (emfAverage==ElectroMotiveForce::arithmetic
                || emfAverage==ElectroMotiveForce::uct0) {
        K_StoreEMF<DIR>(i,j,k,st,sb,Flux,Et,Eb);
      } else if (emfAverage==ElectroMotiveForce::uct_contact) {
        K_StoreContact<DIR>(i,j,k,st,sb,Flux,Et,Eb,SV);
      } else if (emfAverage==ElectroMotiveForce::uct_hll) {
        K_StoreHLL<DIR>(i,j,k,st,sb,sl,sr,vL,vR,Et,Eb,aL,aR,dL,dR);
      }
      /* else if (emfAverage==ElectroMotiveForce::uct_hlld) {
        // We do not have the Alfven speed in the HLL solver
        K_StoreHLLD<DIR>(i,j,k,st,sb,c2Iso,SLb,SRb,vL,vR,uL,uR,Et,Eb,aL,aR,dL,dR);
      }*/
  });

  idfx::popRegion();
}

#endif // HYDRO_MHDSOLVERS_HLLMHD_HPP_
