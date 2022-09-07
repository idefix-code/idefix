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
#include "shockFlattening.hpp"

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
  SlopeLimiter(IdefixArray4D<real> &Vc, IdefixArray1D<real> &dx, ShockFlattening &sf):
        Vc(Vc), dx(dx), flags(sf.flagArray), shockFlattening(sf.isActive) {}

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

  template <typename T>
  KOKKOS_FORCEINLINE_FUNCTION int sign(T val) const {
    return (T(0) < val) - (val < T(0));
}

  KOKKOS_FORCEINLINE_FUNCTION void limitPPMFaceValues(const real vm1, const real v0, const real vp1,
                                                      const real vp2, real &vph) const {
    // if local extremum, then use limited curvature estimate
    if( (vp1-vph)*(vph-v0) < 0.0) {
      // Collela, eqns. 85
      const real deltaL = (vm1-2*v0+vp1);
      const real deltaC = 3*(v0-2*vph+vp1);
      const real deltaR = (v0-2*vp1+vp2);
      // Compute limited curvature estimate
      real delta = 0.0;

      if(sign(deltaL) == sign(deltaC) && sign(deltaR) == sign(deltaC)) {
        const real C = 1.25;
        delta = C * FMIN(FABS(deltaL), FABS(deltaR));
        delta = sign(deltaC) * FMIN(delta, FABS(deltaC));
      }
      vph = 0.5*(v0+vp1) - delta / 6.0;
    }
  }

  // This implementation follows the PPM4 scheme of Peterson & Hammet (PH13)
  // SIAM J. Sci Comput (2013)

  KOKKOS_FORCEINLINE_FUNCTION void getPPMStates(const real vm2, const real vm1, const real v0,
                                                const real vp1, const real vp2, real &vl, real &vr)
                                                const {
    const int n = 2;

    vr = 7.0/12.0*(v0+vp1) - 1.0/12.0*(vm1+vp2);
    vl = 7.0/12.0*(vm1+v0) - 1.0/12.0*(vm2+vp1);

    limitPPMFaceValues(vm2,vm1,v0,vp1,vl);
    limitPPMFaceValues(vm1,v0,vp1,vp2,vr);

    real d2qf = 6.0*(vl + vr - 2.0*v0);
    real d2qc0 = vm1 + vp1 - 2.0*v0;
    real d2qcp1 = v0 + vp2 - 2.0*vp1;
    real d2qcm1 = vm2 + v0 - 2.0*vm1;

    real d2q = 0.0;
    if(sign(d2qf) == sign(d2qc0) && sign(d2qf) == sign(d2qcp1) && sign(d2qf) == sign(d2qcm1)) {
      // smooth extrememum
      const real C = 1.25;
      d2q = FMIN(FABS(d2qc0),FABS(d2qcp1));
      d2q = C * FMIN(FABS(d2qcm1), d2q);
      d2q = sign(d2qf) * FMIN(FABS(d2qf), d2q);
    }

    real qmax = FMAX(FMAX(FABS(vm1),FABS(v0)),FABS(vp1));
    real rho = 0.0;
    // todo(GL): replace 1e-12 by mixed precision value
    if(FABS(d2qf) > 1e-12*qmax) {
      rho = d2q / d2qf;
    }

    // PH13 3.31
    if( ((vr - v0)*(v0 - vl) <= 0) || (vm1 - v0)*(v0 - vp1) <= 0 ) {
      if(rho <= (1.0 - 1e-12)) {
        vl = v0 - rho * (v0 - vl);
        vr = v0 + rho * (vr - v0);
      }
    } else {
      // PH13 3.32
      if(FABS(vr-v0) >= n*FABS(v0-vl)) {
        vr = v0 + n*(v0-vl);
      }
      if(FABS(vl-v0) >= n*FABS(v0-vr)) {
        vl = v0 + n*(v0-vr);
      }
    }
  }

  KOKKOS_FORCEINLINE_FUNCTION void ExtrapolatePrimVar(const int i,
                                                    const int j,
                                                    const int k,
                                                    real vL[], real vR[]) const {
    // 1-- Store the primitive variables on the left, right, and averaged states
    constexpr int ioffset = (dir==IDIR ? 1 : 0);
    constexpr int joffset = (dir==JDIR ? 1 : 0);
    constexpr int koffset = (dir==KDIR ? 1 : 0);

    for(int nv = 0 ; nv < nvmax ; nv++) {
      if constexpr(order == 1) {
        vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset);
        vR[nv] = Vc(nv,k,j,i);
      } else if constexpr(order == 2) {
        real dvm = Vc(nv,k-koffset,j-joffset,i-ioffset)
                  -Vc(nv,k-2*koffset,j-2*joffset,i-2*ioffset);
        real dvp = Vc(nv,k,j,i)-Vc(nv,k-koffset,j-joffset,i-ioffset);

        real dv;
        if(shockFlattening) {
          if(flags(k,j,i) == FlagShock::Shock) {
            // Force slope limiter to minmod
            dv = MinModLim(dvp,dvm);
          } else {
            dv = PLMLim(dvp,dvm);
          }
        } else { // No shock flattening
          dv = PLMLim(dvp,dvm);
        }

        vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset) + HALF_F*dv;

        dvm = dvp;
        dvp = Vc(nv,k+koffset,j+joffset,i+ioffset) - Vc(nv,k,j,i);

        if(shockFlattening) {
          if(flags(k,j,i) == FlagShock::Shock) {
            dv = MinModLim(dvp,dvm);
          } else {
            dv = PLMLim(dvp,dvm);
          }
        } else { // No shock flattening
          dv = PLMLim(dvp,dvm);
        }

        vR[nv] = Vc(nv,k,j,i) - HALF_F*dv;
      } else if constexpr(order == 3) {
          // 1D index along the chosen direction
          const int index = ioffset*i + joffset*j + koffset*k;
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
          real vm2 = Vc(nv,k-3*koffset,j-3*joffset,i-3*ioffset);;
          real vm1 = Vc(nv,k-2*koffset,j-2*joffset,i-2*ioffset);
          real v0 = Vc(nv,k-koffset,j-joffset,i-ioffset);
          real vp1 = Vc(nv,k,j,i);
          real vp2 = Vc(nv,k+koffset,j+joffset,i+ioffset);

          real vr,vl;
          getPPMStates(vm2, vm1, v0, vp1, vp2, vl, vr);
          // vL= left side of current interface (i-1/2)= right side of cell i-1

          // Check positivity
          if(nv==RHO) {
            // If face element is negative, revert to vanleer
            if(vr <= 0.0) {
              real dv = PLMLim(vp1-v0,v0-vm1);
              vr = v0+dv;
            }
          }
          #if HAVE_ENERGY
            if(nv==PRS) {
              // If face element is negative, revert to vanleer
              if(vr <= 0.0) {
                real dv = PLMLim(vp1-v0,v0-vm1);
                vr = v0+dv;
              }
            }
          #endif

          vL[nv] = vr;
          // Reconstruction in cell i

          vm2 = vm1;
          vm1 = v0;
          v0 = vp1;
          vp1 = vp2;
          vp2 = Vc(nv,k+2*koffset,j+2*joffset,i+2*ioffset);

          getPPMStates(vm2, vm1, v0, vp1, vp2, vl, vr);

          // Check positivity
          if(nv==RHO) {
            // If face element is negative, revert to vanleer
            if(vl <= 0.0) {
              real dv = PLMLim(vp1-v0,v0-vm1);
              vl = v0-dv;
            }
          }
          #if HAVE_ENERGY
            if(nv==PRS) {
              // If face element is negative, revert to vanleer
              if(vl <= 0.0) {
                real dv = PLMLim(vp1-v0,v0-vm1);
                vl = v0-dv;
              }
            }
          #endif

          vR[nv] = vl;
      }
    }
  }

  IdefixArray4D<real> Vc;
  IdefixArray1D<real> dx;
  IdefixArray3D<FlagShock> flags;
  bool shockFlattening{false};
};


#endif // HYDRO_SLOPELIMITER_HPP_
