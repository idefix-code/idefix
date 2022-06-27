// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************
#ifndef HYDRO_SLOPELIMITER_HPP_
#define HYDRO_SLOPELIMITER_HPP_

#include "hydro.hpp"
#include "dataBlock.hpp"

// Build a left and right extrapolation of the primitive variables along direction dir

// These functions extrapolate the cell prim vars to the faces. Definitions are as followed
//
// |       cell i-1               interface i          cell i
// |-----------------------------------|------------------------------------||
//          Vc(i-1)           PrimL(i)  PrimR(i)       Vc(i)
template<const int dir,
         const int nvmax,
         const Limiter limiter = Limiter::VanLeer,
         const int order = ORDER>
class SlopeLimiter {
 public:
  SlopeLimiter(IdefixArray4D<real> &Vc, IdefixArray1D<real> &dx): Vc(Vc), dx(dx) {}

  KOKKOS_FORCEINLINE_FUNCTION real MinModLim(const real dvp, const real dvm) const {
    real dq= 0.0;
    // MinMod
    if(dvp*dvm >0.0) {
      real dq = ( fabs(dvp) < fabs(dvm) ? dvp : dvm);
    }
    return(dq);
  }

  KOKKOS_FORCEINLINE_FUNCTION real LimO3Lim(const real dvp, const real dvm, const real dx) const {
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

  KOKKOS_FORCEINLINE_FUNCTION real VanLeerLim(const real dvp, const real dvm) const {
    real dq= 0.0;
    dq = (dvp*dvm > ZERO_F ? TWO_F*dvp*dvm/(dvp + dvm) : ZERO_F);
    return(dq);
  }

  KOKKOS_FORCEINLINE_FUNCTION real McLim(const real dvp, const real dvm) const {
    real dq = 0;
    if(dvp*dvm >0.0) {
      real dqc = 0.5*(dvp+dvm);
      real d2q = 2.0*( fabs(dvp) < fabs(dvm) ? dvp : dvm);
      dq= fabs(d2q) < fabs(dqc) ? d2q : dqc;
    }
    return(dq);
  }

  KOKKOS_FORCEINLINE_FUNCTION real PLMLim(const real dvp, const real dvm) const {
    if constexpr(limiter == Limiter::VanLeer) return(VanLeerLim(dvp,dvm));
    if constexpr(limiter == Limiter::McLim) return(McLim(dvp,dvm));
    if constexpr(limiter == Limiter::MinMod) return(MinModLim(dvp,dvm));
  }

  KOKKOS_FORCEINLINE_FUNCTION void ExtrapolatePrimVar(const int i,
                                                    const int j,
                                                    const int k,
                                                    real vL[], real vR[]) const {
    // 1-- Store the primitive variables on the left, right, and averaged states
    constexpr int ioffset = (dir==IDIR ? 1 : 0);
    constexpr int joffset = (dir==JDIR ? 1 : 0);
    constexpr int koffset = (dir==KDIR ? 1 : 0);

    // 1D index along the chosen direction
    const int index = ioffset*i + joffset*j + koffset*k;
    for(int nv = 0 ; nv < nvmax ; nv++) {
      if constexpr(order == 1) {
        vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset);
        vR[nv] = Vc(nv,k,j,i);
      } else if constexpr(order == 2) {
        real dvm = Vc(nv,k-koffset,j-joffset,i-ioffset)
                  -Vc(nv,k-2*koffset,j-2*joffset,i-2*ioffset);
        real dvp = Vc(nv,k,j,i)-Vc(nv,k-koffset,j-joffset,i-ioffset);
        real dv = PLMLim(dvp,dvm);

        vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset) + HALF_F*dv;

        dvm = dvp;
        dvp = Vc(nv,k+koffset,j+joffset,i+ioffset) - Vc(nv,k,j,i);
        dv = PLMLim(dvp,dvm);

        vR[nv] = Vc(nv,k,j,i) - HALF_F*dv;
      } else if constexpr(order == 3) {
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
      } else if constexpr(order == 4) {
          // Reconstruction in cell i-1
          real dvm2 = Vc(nv,k-2*koffset,j-2*joffset,i-2*ioffset)
                    -Vc(nv,k-3*koffset,j-3*joffset,i-3*ioffset);
          real dvm1 = Vc(nv,k-koffset,j-joffset,i-ioffset)
                    -Vc(nv,k-2*koffset,j-2*joffset,i-2*ioffset);
          real dvp1 = Vc(nv,k,j,i) - Vc(nv,k-koffset,j-joffset,i-ioffset);
          real dvp2 = Vc(nv,k+koffset,j+joffset,i+ioffset) - Vc(nv,k,j,i);

          real dvlm = McLim(dvm1,dvm2);
          real dvl0 = McLim(dvp1,dvm1);
          real dvlp = McLim(dvp2,dvp1);

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
          dvlp = McLim(dvp2,dvp1);

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
    }
  }

  IdefixArray4D<real> &Vc;
  IdefixArray1D<real> &dx;
};


#endif // HYDRO_SLOPELIMITER_HPP_
