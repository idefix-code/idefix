// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_RIEMANNSOLVER_RADIATIONSOLVERS_LIM_FLUXRAD_HPP_
#define FLUID_RIEMANNSOLVER_RADIATIONSOLVERS_LIM_FLUXRAD_HPP_

//
#include <string>
#include <tuple>
#include "idefix.hpp"
#include "input.hpp"


KOKKOS_INLINE_FUNCTION void K_LimitRadFluxSimple(real  V[]) {
real Fnorm = std::sqrt(EXPAND(V[FR1]*V[FR1] , + V[FR2]*V[FR2], + V[FR3]*V[FR3]));

if (Fnorm > V[ER]) {
    EXPAND( V[FR1] *= (Fnorm <= 1.e-50 ? V[ER]/1.e-50 : (V[ER]/Fnorm)); ,
            V[FR2] *= (Fnorm <= 1.e-50 ? V[ER]/1.e-50 : (V[ER]/Fnorm)); ,
            V[FR3] *= (Fnorm <= 1.e-50 ? V[ER]/1.e-50 : (V[ER]/Fnorm)); )
}

return;
}

KOKKOS_INLINE_FUNCTION void K_LimitRadFluxFlatten(
    real  vL[], real vR[], real Vc[], real VOff[]) {
real FnormL = std::sqrt(EXPAND(vL[FR1]*vL[FR1] , + vL[FR2]*vL[FR2], + vL[FR3]*vL[FR3]));
real FnormR = std::sqrt(EXPAND(vR[FR1]*vR[FR1] , + vR[FR2]*vR[FR2], + vR[FR3]*vR[FR3]));

real Fnorm = std::sqrt(EXPAND(Vc[FR1]*Vc[FR1] , + Vc[FR2]*Vc[FR2], + Vc[FR3]*Vc[FR3]));
if (Fnorm > Vc[ER]) {
    //printf("Limit Flux Fnorm/Er=%e \n",Fnorm/V[ER]);
    EXPAND( Vc[FR1] *= (Fnorm <= 1.e-50 ? Vc[ER]/1.e-50 : Vc[ER]/Fnorm); ,
            Vc[FR2] *= (Fnorm <= 1.e-50 ? Vc[ER]/1.e-50 : Vc[ER]/Fnorm); ,
            Vc[FR3] *= (Fnorm <= 1.e-50 ? Vc[ER]/1.e-50 : Vc[ER]/Fnorm); )
}

real FnormOff = std::sqrt( EXPAND( VOff[FR1]*VOff[FR1] ,
                                      + VOff[FR2]*VOff[FR2],
                                      + VOff[FR3]*VOff[FR3]) );
if (FnormOff > VOff[ER]) {
    //printf("Limit Flux Fnorm/Er=%e \n",Fnorm/V[ER]);
    EXPAND( VOff[FR1] *= (FnormOff <= 1.e-50 ? VOff[ER]/1.e-50 : VOff[ER]/FnormOff); ,
            VOff[FR2] *= (FnormOff <= 1.e-50 ? VOff[ER]/1.e-50 : VOff[ER]/FnormOff); ,
            VOff[FR3] *= (FnormOff <= 1.e-50 ? VOff[ER]/1.e-50 : VOff[ER]/FnormOff); )
}

if ((FnormL > vL[ER]) || (FnormR > vR[ER])) {
    //printf("Limit Flux Fnorm/Er=%e \n",Fnorm/V[ER]);
    EXPAND( vL[ER] = VOff[ER];
            vR[ER] = Vc[ER];
            vL[FR1] = VOff[FR1];
            vR[FR1] = Vc[FR1]; ,
            vL[FR2] = VOff[FR2];
            vR[FR2] = Vc[FR2]; ,
            vL[FR3] = VOff[FR3];
            vR[FR3] = Vc[FR3]; )
}

return;
}

KOKKOS_INLINE_FUNCTION void K_LimitRadFluxReconstruct(
    real  vL[], real vR[], real Vc[], real VOff[]) {
real FnormL = std::sqrt(EXPAND(vL[FR1]*vL[FR1] , + vL[FR2]*vL[FR2], + vL[FR3]*vL[FR3]));
real FnormR = std::sqrt(EXPAND(vR[FR1]*vR[FR1] , + vR[FR2]*vR[FR2], + vR[FR3]*vR[FR3]));
real FnormOff = std::sqrt( EXPAND( VOff[FR1]*VOff[FR1] ,
                                   + VOff[FR2]*VOff[FR2],
                                   + VOff[FR3]*VOff[FR3]));
real Fnorm = std::sqrt(EXPAND(Vc[FR1]*Vc[FR1] , + Vc[FR2]*Vc[FR2], + Vc[FR3]*Vc[FR3]));

real ratiocellOff = FnormOff/VOff[ER];
real ratiocell = Fnorm/Vc[ER];
real ratioL = FnormL/vL[ER];
real ratioR = FnormR/vR[ER];

if ((FnormL > vL[ER]) || (ratiocellOff > ratioL)) {
    EXPAND( vL[FR1] = VOff[FR1]*vL[ER]/VOff[ER]; ,
            vL[FR2] = VOff[FR2]*vL[ER]/VOff[ER]; ,
            vL[FR3] = VOff[FR3]*vL[ER]/VOff[ER]; )
}
if ((FnormR > vR[ER])|| (ratiocell > ratioR)) {
    EXPAND( vR[FR1] = Vc[FR1]*vR[ER]/Vc[ER]; ,
            vR[FR2] = Vc[FR2]*vR[ER]/Vc[ER]; ,
            vR[FR3] = Vc[FR3]*vR[ER]/Vc[ER]; )
}


return;
}

template <typename Phys>
KOKKOS_INLINE_FUNCTION void K_LimitRadFlux(
    real  vL[], real vR[], real Vc[], real VOff[],
    typename RiemannSolver<Phys>::Type_Radlimiter radLimiter) {

  if (radLimiter == RiemannSolver<Phys>::Type_Radlimiter::simple) {
      K_LimitRadFluxSimple(vL);
      K_LimitRadFluxSimple(vR);
  } else if (radLimiter == RiemannSolver<Phys>::Type_Radlimiter::flatten) {
      K_LimitRadFluxFlatten(vL,vR,Vc,VOff);
  } else if (radLimiter == RiemannSolver<Phys>::Type_Radlimiter::fpreserving) {
      K_LimitRadFluxReconstruct(vL,vR,Vc,VOff);
  }

return;
}


#endif //FLUID_RIEMANNSOLVER_RADIATIONSOLVERS_LIM_FLUXRAD_HPP_
