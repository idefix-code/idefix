// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************
#ifndef FLUID_RIEMANNSOLVER_SLOPELIMITER_HPP_
#define FLUID_RIEMANNSOLVER_SLOPELIMITER_HPP_

// The default PLM Limiter template type
enum class PLMLimiter {VanLeer, MinMod, McLim};

template<const PLMLimiter limiter = PLMLimiter::VanLeer>
class SlopeLimiter {
 public:
  KOKKOS_FORCEINLINE_FUNCTION static real MinModLim(const real dvp, const real dvm) {
    real dq= 0.0;
    // MinMod
    if(dvp*dvm >0.0) {
      real dq = ( fabs(dvp) < fabs(dvm) ? dvp : dvm);
    }
    return(dq);
  }

  KOKKOS_FORCEINLINE_FUNCTION real static LimO3Lim(const real dvp, const real dvm, const real dx) {
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

  KOKKOS_FORCEINLINE_FUNCTION static real VanLeerLim(const real dvp, const real dvm) {
    real dq = (dvp*dvm > ZERO_F ? TWO_F*dvp*dvm/(dvp + dvm) : ZERO_F);
    return(dq);
  }

  // Generalize vanleer for non-homogeneous grids
  KOKKOS_FORCEINLINE_FUNCTION static real VanLeerLim(const real dvp, const real dvm,
                                              const real cp , const real cm) {
    real dq = (dvp*dvm > 0.0 ? dvp*dvm*(cp*dvm + cm*dvp)
                       /(dvp*dvp + dvm*dvm + (cp + cm - 2.0)*dvp*dvm) : 0.0);
    return dq;
  }

  KOKKOS_FORCEINLINE_FUNCTION static real McLim(const real dvp, const real dvm) {
    real dq = 0;
    if(dvp*dvm >0.0) {
      real dqc = 0.5*(dvp+dvm);
      real d2q = 2.0*( fabs(dvp) < fabs(dvm) ? dvp : dvm);
      dq= fabs(d2q) < fabs(dqc) ? d2q : dqc;
    }
    return(dq);
  }

  // Generalized McLimiter for non-homogeneous grid
  KOKKOS_FORCEINLINE_FUNCTION static real McLim(const real dvp, const real dvm,
                                          const real cp , const real cm) {
    real dq = 0;
    if(dvp*dvm >0.0) {
      real dqc = 0.5*(dvp+dvm);
      real d2q =  fabs(dvp*cp) < fabs(dvm*cm) ? dvp*cp : dvm*cm;
      dq= fabs(d2q) < fabs(dqc) ? d2q : dqc;
    }
    return(dq);
  }


  KOKKOS_FORCEINLINE_FUNCTION static real PLMLim(const real dvp, const real dvm) {
    if constexpr(limiter == PLMLimiter::VanLeer) return(VanLeerLim(dvp,dvm));
    if constexpr(limiter == PLMLimiter::McLim) return(McLim(dvp,dvm));
    if constexpr(limiter == PLMLimiter::MinMod) return(MinModLim(dvp,dvm));
  }

  // Overlad of PLM limiter for irregular grids
  KOKKOS_FORCEINLINE_FUNCTION static real PLMLim(const real dvp, const real dvm,
                                          const real cp, const real cm) {
    if constexpr(limiter == PLMLimiter::VanLeer) return(VanLeerLim(dvp,dvm,cp,cm));
    if constexpr(limiter == PLMLimiter::McLim) return(McLim(dvp,dvm,cp,cm));
    if constexpr(limiter == PLMLimiter::MinMod) return(MinModLim(dvp,dvm));
  }


  template <typename T>
  KOKKOS_FORCEINLINE_FUNCTION static int sign(T val) {
    return (T(0) < val) - (val < T(0));
  }

  // PPM limiter, inspired from
  // PH13: Peterson, J. L. & Hammett, G. W. Positivity Preservation and Advection Algorithms
  //       with Applications to Edge Plasma Turbulence. SIAM J. Sci. Comput. 35, B576–B605 (2013).
  // CD11: Colella, P., Dorr, M. R., Hittinger, J. A. F. & Martin, D. F. High-order,
  //       finite-volume methods in mapped coordinates. Journal of Computational Physics 230,
  //       2952–2976 (2011).
  // CS08: Colella, P. & Sekora, M. D. A limiter for PPM that preserves accuracy at smooth extrema.
  //       Journal of Computational Physics 227, 7069–7076 (2008).
  // FS18: Felker, K. G. & Stone, J. M. A fourth-order accurate finite volume method for ideal MHD
  //       via upwind constrained transport. Journal of Computational Physics 375, 1365–1400 (2018).

  KOKKOS_FORCEINLINE_FUNCTION static void limitPPMFaceValues(const real vm1, const real v0,
                                                      const real vp1, const real vp2, real &vph) {
    // if local extremum, then use limited curvature estimate
    if( (vp1-vph)*(vph-v0) < 0.0) {
      // CD11, eqns. 85
      const real deltaL = (vm1-2*v0+vp1);
      const real deltaC = 3*(v0-2*vph+vp1);
      const real deltaR = (v0-2*vp1+vp2);
      // Compute limited curvature estimate
      real delta = 0.0;

      // CS08 eq. 18 with corrections from FS18 section. 2.2.2
      if(sign(deltaL) == sign(deltaC) && sign(deltaR) == sign(deltaC)) {
        const real C = 1.25;
        delta = C * FMIN(FABS(deltaL), FABS(deltaR));
        delta = sign(deltaC) * FMIN(delta, FABS(deltaC));
      }
      vph = 0.5*(v0+vp1) - delta / 6.0; // CD11 eq. 88 (correction of CS08, eq. 19)
    }
  }

  KOKKOS_FORCEINLINE_FUNCTION static void getPPMStates(const real vm2, const real vm1,
                                                      const real v0, const real vp1, const real vp2,
                                                      real &vl, real &vr) {
    const int n = 2;

    // 1: unlimited left and right interpolant (PH13 3.26-3.27)
    vr = 7.0/12.0*(v0+vp1) - 1.0/12.0*(vm1+vp2);
    vl = 7.0/12.0*(vm1+v0) - 1.0/12.0*(vm2+vp1);

    // 2: limit interpolated face values (CD11 4.3.1)
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
};

#endif // FLUID_RIEMANNSOLVER_SLOPELIMITER_HPP_
