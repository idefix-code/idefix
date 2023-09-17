// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************
#ifndef FLUID_BRAGINSKII_BRAGLIMITER_HPP_
#define FLUID_BRAGINSKII_BRAGLIMITER_HPP_

#include "idefix.hpp"

template<const Limiter limiter = Limiter::VanLeer>
class BragLimiter {
 public:
  BragLimiter()=default;

  KOKKOS_FORCEINLINE_FUNCTION static real MinModLim(const real dvp, const real dvm) {
    real dq= 0.0;
    // MinMod
    if(dvp*dvm >0.0) {
      real dq = ( fabs(dvp) < fabs(dvm) ? dvp : dvm);
    }
    return(dq);
  }

  KOKKOS_FORCEINLINE_FUNCTION static real VanLeerLim(const real dvp, const real dvm) {
    real dq= 0.0;
    dq = (dvp*dvm > ZERO_F ? TWO_F*dvp*dvm/(dvp + dvm) : ZERO_F);
    return(dq);
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

  KOKKOS_FORCEINLINE_FUNCTION static real Lim(const real dvp, const real dvm) {
    if constexpr(limiter == Limiter::VanLeer) return(VanLeerLim(dvp,dvm));
    if constexpr(limiter == Limiter::McLim) return(McLim(dvp,dvm));
    if constexpr(limiter == Limiter::MinMod) return(MinModLim(dvp,dvm));
  }
};

#endif //FLUID_BRAGINSKII_BRAGLIMITER_HPP_
