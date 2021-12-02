// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************
#ifndef HYDRO_EXTRAPOLATEPRIMVAR_HPP_
#define HYDRO_EXTRAPOLATEPRIMVAR_HPP_

#include "hydro.hpp"
#include "dataBlock.hpp"

// LimO3 Slope limiter
KOKKOS_FORCEINLINE_FUNCTION real LimO3Lim(real dvp, real dvm, real dx) {
  real r = 0.1;
  real a,b,c,q, th, lim;
  real eta, psi, eps = 1.e-12;

  th  = dvm/(dvp + 1.e-16);

  q = (2.0 + th)/3.0;

  a = FMIN(1.5,2.0*th);
  a = FMIN(q,a);
  b = FMAX(-0.5*th,a);
  c = FMIN(q,b);
  psi = FMAX(0.0,c);

  eta = r*dx;
  eta = (dvm*dvm + dvp*dvp)/(eta*eta);
  if ( eta <= 1.0 - eps) {
    lim = q;
  } else if (eta >= 1.0 + eps) {
    lim = psi;
  } else {
    psi =   (1.0 - (eta - 1.0)/eps)*q
          + (1.0 + (eta - 1.0)/eps)*psi;
    lim = 0.5*psi;
  }
  return (lim);
}

KOKKOS_FORCEINLINE_FUNCTION real PPMLim(real dvp, real dvm) {
  if(dvp*dvm >0.0) {
    real dqc = 0.5*(dvp+dvm);
    real d2q = 2.0*( fabs(dvp) < fabs(dvm) ? dvp : dvm);
    return( fabs(d2q) < fabs(dqc) ? d2q : dqc);
  }
  return(ZERO_F);
}

// Build a left and right extrapolation of the primitive variables along direction dir

// These functions extrapolate the cell prim vars to the faces. Definitions are as followed
//
// |       cell i-1               interface i          cell i
// |-----------------------------------|------------------------------------||
//          Vc(i-1)           PrimL(i)  PrimR(i)       Vc(i)

template<const int DIR>
KOKKOS_FORCEINLINE_FUNCTION void K_ExtrapolatePrimVar
      (const int i, const int j, const int k, const IdefixArray4D<real> &Vc,
      const IdefixArray4D<real> &Vs, const IdefixArray1D<real> &dx, real vL[], real vR[]) {
  // 1-- Store the primitive variables on the left, right, and averaged states
  const int ioffset = (DIR==IDIR ? 1 : 0);
  const int joffset = (DIR==JDIR ? 1 : 0);
  const int koffset = (DIR==KDIR ? 1 : 0);
  const int BXn = BX1 + DIR;

  // 1D index along the chosen direction
  const int index = ioffset*i + joffset*j + koffset*k;

  #pragma unroll
  for(int nv = 0 ; nv < NVAR; nv++) {
    #if ORDER == 1
      if(nv==BXn) {
        vR[nv] = vL[nv] = Vs(DIR,k,j,i);
      } else {
        vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset);
        vR[nv] = Vc(nv,k,j,i);
      }
    #elif ORDER == 2
      if(nv==BXn) {
        vR[nv] = vL[nv] = Vs(DIR,k,j,i);
      } else {
        real dvm = Vc(nv,k-koffset,j-joffset,i-ioffset)
                  -Vc(nv,k-2*koffset,j-2*joffset,i-2*ioffset);
        real dvp = Vc(nv,k,j,i)-Vc(nv,k-koffset,j-joffset,i-ioffset);
        // Van Leer limiter
        real dv = (dvp*dvm > ZERO_F ? TWO_F*dvp*dvm/(dvp + dvm) : ZERO_F);

        vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset) + HALF_F*dv;

        dvm = dvp;
        dvp = Vc(nv,k+koffset,j+joffset,i+ioffset) - Vc(nv,k,j,i);

        // Van Leer limiter
        dv = (dvp*dvm > ZERO_F ? TWO_F*dvp*dvm/(dvp + dvm) : ZERO_F);

        vR[nv] = Vc(nv,k,j,i) - HALF_F*dv;
      }
    #elif ORDER == 3
      if(nv==BXn) {
        vR[nv] = vL[nv] = Vs(DIR,k,j,i);
      } else {
        real dvm = Vc(nv,k-koffset,j-joffset,i-ioffset)
                  -Vc(nv,k-2*koffset,j-2*joffset,i-2*ioffset);
        real dvp = Vc(nv,k,j,i)-Vc(nv,k-koffset,j-joffset,i-ioffset);

        // Limo3 limiter
        real dv = dvp * LimO3Lim(dvp, dvm, dx(index-1));
        vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset) + HALF_F*dv;

        // Check positivity
        if(nv==RHO) {
          // If face element is negative, revert to vanleer
          if(vL[nv] <= 0.0) {
            dv = (dvp*dvm > ZERO_F ? TWO_F*dvp*dvm/(dvp + dvm) : ZERO_F);
            vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset) + HALF_F*dv;
          }
        }
        #if HAVE_ENERGY
          if(nv==PRS) {
            // If face element is negative, revert to vanleer
            if(vL[nv] <= 0.0) {
              dv = (dvp*dvm > ZERO_F ? TWO_F*dvp*dvm/(dvp + dvm) : ZERO_F);
              vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset) + HALF_F*dv;
            }
          }
        #endif

        dvm = dvp;
        dvp = Vc(nv,k+koffset,j+joffset,i+ioffset) - Vc(nv,k,j,i);

        // Limo3 limiter
        dv = dvm * LimO3Lim(dvm, dvp, dx(index));
        vR[nv] = Vc(nv,k,j,i) - HALF_F*dv;

        // Check positivity
        if(nv==RHO) {
          // If face element is negative, revert to vanleer
          if(vR[nv] <= 0.0) {
            dv = (dvp*dvm > ZERO_F ? TWO_F*dvp*dvm/(dvp + dvm) : ZERO_F);
            vR[nv] = Vc(nv,k,j,i) - HALF_F*dv;
          }
        }
        #if HAVE_ENERGY
          if(nv==PRS) {
            // If face element is negative, revert to vanleer
            if(vR[nv] <= 0.0) {
              dv = (dvp*dvm > ZERO_F ? TWO_F*dvp*dvm/(dvp + dvm) : ZERO_F);
              vR[nv] = Vc(nv,k,j,i) - HALF_F*dv;
            }
          }
        #endif
      }
    #elif ORDER == 4
    if(nv==BXn) {
        vR[nv] = vL[nv] = Vs(DIR,k,j,i);
      } else {
        // Reconstruction in cell i-1
        real dvm2 = Vc(nv,k-2*koffset,j-2*joffset,i-2*ioffset)
                  -Vc(nv,k-3*koffset,j-3*joffset,i-3*ioffset);
        real dvm1 = Vc(nv,k-koffset,j-joffset,i-ioffset)
                  -Vc(nv,k-2*koffset,j-2*joffset,i-2*ioffset);
        real dvp1 = Vc(nv,k,j,i) - Vc(nv,k-koffset,j-joffset,i-ioffset);
        real dvp2 = Vc(nv,k+koffset,j+joffset,i+ioffset) - Vc(nv,k,j,i);

        real dvlm = PPMLim(dvm1,dvm2);
        real dvl0 = PPMLim(dvp1,dvm1);
        real dvlp = PPMLim(dvp2,dvp1);

        real dvp = 0.5 * dvp1 - (dvlp - dvl0) / 6.0;
        real dvm = -0.5 * dvm1 - (dvl0 - dvlm) / 6.0;

        real vc = Vc(nv,k-koffset,j-joffset,i-ioffset);

        real vl = vc + dvm;
        real vr = vc + dvp;

        if(dvp*dvm>0.0) {
          vl = vr = vc;
        } else {
          if( (vr - vl) * (vc-0.5*(vl+vr)) > (vr-vl)*(vr-vl) / 6.0) {
            vl = 3.0*vc - 2.0*vr;
          }
          if( (vr - vl) * (vc-0.5*(vl+vr)) < - (vr-vl)*(vr-vl) / 6.0) {
            vr = 3.0*vc - 2.0*vl;
          }

          //if(FABS(dvp) >= 2.0*FABS(dvm)) dvp = -2.0*dvm;
          // Skip because it's not used
          //if(FABS(dqm) >= 2.0*FABS(dqp)) dqm = -2.0*dqp;
        }

        // vL= left side of current interface (i-1/2)= right side of cell i-1
        vL[nv] = vr;

        // Reconstruction in cell i

        dvm2 = dvm1;
        dvm1 = dvp1;
        dvp1 = dvp2;
        dvp2 = Vc(nv,k+2*koffset,j+2*joffset,i+2*ioffset)
                - Vc(nv,k+koffset,j+joffset,i+ioffset);

        dvlm = dvl0;
        dvl0 = dvlp;
        dvlp = PPMLim(dvp2,dvp1);

        dvp = 0.5 * dvp1 - (dvlp - dvl0) / 6.0;
        dvm = -0.5 * dvm1 - (dvl0 - dvlm) / 6.0;

        vc = Vc(nv,k,j,i);

        vl = vc + dvm;
        vr = vc + dvp;

        if(dvp*dvm>0.0) {
          vl = vr = vc;
        } else {
          if( (vr - vl) * (vc-0.5*(vl+vr)) > (vr-vl)*(vr-vl) / 6.0) {
            vl = 3.0*vc - 2.0*vr;
          }
          if( (vr - vl) * (vc-0.5*(vl+vr)) < - (vr-vl)*(vr-vl) / 6.0) {
            vr = 3.0*vc - 2.0*vl;
          }

          //if(FABS(dvp) >= 2.0*FABS(dvm)) dvp = -2.0*dvm;
          // Skip because it's not used
          //if(FABS(dqm) >= 2.0*FABS(dqp)) dqm = -2.0*dqp;
        }
        // vR= right side of current interface (i-1/2)= left side of cell i
        vR[nv] = vl;
      }
    #endif
  }
}
#endif // HYDRO_EXTRAPOLATEPRIMVAR_HPP_
