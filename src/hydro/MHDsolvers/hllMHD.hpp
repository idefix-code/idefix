// ***********************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************************

#ifndef HYDRO_MHDSOLVERS_HLLMHD_HPP_
#define HYDRO_MHDSOLVERS_HLLMHD_HPP_

#include "../idefix.hpp"
#include "solversMHD.hpp"

// Compute Riemann fluxes from states using HLL solver
template<const int DIR, const int Xn, const int Xt, const int Xb,
         const int BXn, const int BXt, const int BXb>
void HllMHD(DataBlock & data, real gamma, real C2Iso, ParabolicType haveHallin, real xH) {
  idfx::pushRegion("HLL_MHD");
  
  int ioffset,joffset,koffset;
  int iextend, jextend,kextend;
  ioffset=joffset=koffset=0;
  // extension in perp to the direction of integration, as required by CT.
  iextend=jextend=kextend=0;

  IdefixArray4D<real> PrimL = data.PrimL;
  IdefixArray4D<real> PrimR = data.PrimR;
  IdefixArray4D<real> Flux = data.FluxRiemann;
  IdefixArray3D<real> cMax = data.cMax;

  ParabolicType haveHall = haveHallin;
  IdefixArray4D<real> J = data.J;
  IdefixArray3D<real> xHallArr = data.xHall;
  IdefixArray1D<real> dx = data.dx[DIR];
  IdefixArray1D<real> dx2 = data.dx[JDIR];
  IdefixArray1D<real> x1 = data.x[IDIR];
  IdefixArray1D<real> rt = data.rt;
  IdefixArray1D<real> dmu = data.dmu;

  // References to required emf components
  IdefixArray3D<real> Eb;
  IdefixArray3D<real> Et;
  
  IdefixArray3D<int> SV;

  real xHConstant = xH;
  real gamma_m1=gamma-ONE_F;

  // Define normal, tangent and bi-tanget indices
  // st and sb will be useful only when Hall is included
  D_EXPAND( real st;  ,
                      ,
            real sb;  )

  switch(DIR) {
    case(IDIR):
      ioffset = 1;
      D_EXPAND(               ,
                jextend = 1;  ,
                kextend = 1;  )

      Et = data.emf.ezi;
      Eb = data.emf.eyi;
      SV = data.emf.svx;

      D_EXPAND( st = -1.0;  ,
                            ,
                sb = +1.0;  )
      break;
    case(JDIR):
      joffset=1;
      D_EXPAND( iextend = 1;  ,
                              ,
                kextend = 1;  )

      Et = data.emf.ezj;
      Eb = data.emf.exj;
      SV = data.emf.svy;

      D_EXPAND( st = +1.0;  ,
                            ,
                sb = -1.0;  )
      break;
    case(KDIR):
      koffset=1;
      D_EXPAND( iextend = 1;  ,
              jextend = 1;    ,
              )
      
      Et = data.emf.eyk;
      Eb = data.emf.exk;
      SV = data.emf.svz;

      D_EXPAND( st = -1.0;  ,
                            ,
                sb = +1.0;  )
      break;
    default:
      IDEFIX_ERROR("Wrong direction");
  }


  idefix_for("CalcRiemannFlux",
             data.beg[KDIR]-kextend,data.end[KDIR]+koffset+kextend,
             data.beg[JDIR]-jextend,data.end[JDIR]+joffset+jextend,
             data.beg[IDIR]-iextend,data.end[IDIR]+ioffset+iextend,
    KOKKOS_LAMBDA (int k, int j, int i) {
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
      real gpr, b1, b2, b3, Btmag2, Bmag2;
      real xH;
#if HAVE_ENERGY
      gpr = gamma*vL[PRS];
#else
      gpr = C2Iso*vL[RHO];
#endif

      // -- get total field
      b1 = b2 = b3 = ZERO_F;
      EXPAND (b1 = vL[BXn];  ,
              b2 = vL[BXt];  ,
              b3 = vL[BXb];)

      Btmag2 = b2*b2 + b3*b3;
      Bmag2  = b1*b1 + Btmag2;

      cL = gpr - Bmag2;
      cL = gpr + Bmag2 + sqrt(cL*cL + 4.0*gpr*Btmag2);
      cL = sqrt(HALF_F*cL/vL[RHO]);
      
#if HAVE_ENERGY
      gpr = gamma*vR[PRS];
#else
      gpr = C2Iso*vR[RHO];
#endif

      // -- get total field
      b1 = b2 = b3 = ZERO_F;
      EXPAND (b1 = vR[BXn];  ,
              b2 = vR[BXt];  ,
              b3 = vR[BXb];)

      Btmag2 = b2*b2 + b3*b3;
      Bmag2  = b1*b1 + Btmag2;

      cR = gpr - Bmag2;
      cR = gpr + Bmag2 + sqrt(cR*cR + 4.0*gpr*Btmag2);
      cR = sqrt(HALF_F*cR/vR[RHO]);
      
      // 4.1 
      real cminL = vL[Xn] - cL;
      real cmaxL = vL[Xn] + cL;
      
      real cminR = vR[Xn] - cR;
      real cmaxR = vR[Xn] + cR;
      
      real SL = FMIN(cminL, cminR);
      real SR = FMAX(cmaxL, cmaxR);
      
      // Signal speeds specific to B (different from the other ones when Hall is enabled)
      real SLb = SL;
      real SRb = SR;
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

        real cw = FABS(xH) * sqrt(Bmag2) / dl;

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
      K_Flux(fluxL, vL, fluxL, C2Iso, Xn, Xt, Xb, BXn, BXt, BXb);
      K_Flux(fluxR, vR, fluxR, C2Iso, Xn, Xt, Xb, BXn, BXt, BXb);

      // 4-- Compute the Hall flux
      if(haveHall) {
        int ip1, jp1, kp1;
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

          fluxL[BX2] += -xH* (  Jx3*uL[BX2]  );
          fluxR[BX2] += -xH* (  Jx3*uR[BX2]  );

          #if COMPONENTS==3
          fluxL[BX1] += -xH* (   - Jx1*uL[BX3] );
          fluxR[BX1] += -xH* (   - Jx1*uR[BX3] );

          fluxL[BX2] += -xH* (   - Jx2*uL[BX3] );
          fluxR[BX2] += -xH* (   - Jx2*uR[BX3] );
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
      if (SL > 0) {
        #pragma unroll
        for (int nv = 0 ; nv < NFLX; nv++) {
          Flux(nv,k,j,i) = fluxL[nv];
        }
      } else if (SR < 0) {
        #pragma unroll
        for (int nv = 0 ; nv < NFLX; nv++) {
          Flux(nv,k,j,i) = fluxR[nv];
        }
      } else {
        #pragma unroll
        for(int nv = 0 ; nv < NFLX; nv++) {
          Flux(nv,k,j,i) = SL*SR*uR[nv] - SL*SR*uL[nv] + SR*fluxL[nv] - SL*fluxR[nv];
          Flux(nv,k,j,i) *= (1.0 / (SR - SL));
        }
      }

      //6-- Compute maximum wave speed for this sweep
      cMax(k,j,i) = cmax;

      // 7-- Store the flux in the emf components
      D_EXPAND( Et(k,j,i) = st*Flux(BXt,k,j,i); ,
                                                ,
                Eb(k,j,i) = sb*Flux(BXb,k,j,i); )
      
#if EMF_AVERAGE == UCT_CONTACT
      int s = 0;
      if (Flux(RHO,k,j,i) >  eps_UCT_CONTACT) s =  1;
      if (Flux(RHO,k,j,i) < -eps_UCT_CONTACT) s = -1;

      SV(k,j,i) = s;
#endif
  });

  idfx::popRegion();
}

#endif // HYDRO_MHDSOLVERS_HLLMHD_HPP_
