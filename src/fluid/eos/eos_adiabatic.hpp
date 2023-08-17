// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_EOS_EOS_ADIABATIC_HPP_
#define FLUID_EOS_EOS_ADIABATIC_HPP_

#include <string>
#include "idefix.hpp"
#include "input.hpp"

// This is the adiabatic (constant gamma) implementation of the equation of state
class EquationOfState {
 public:
  EquationOfState() = default;

  EquationOfState(Input & input, DataBlock *, std::string prefix) {
    this->gamma = input.GetOrSet<real>(prefix,"gamma",0, 5.0/3.0);
  }

  void ShowConfig() {
    idfx::cout << "EquationOfState: ideal with gamma=" << this->gamma << std::endl;
  }

  // First adiabatic exponent. In the ideal EOS, gamma does not depend on the gas state,
  // So we add default values to 0 here so that GetGamma can be called without any argument
  KOKKOS_INLINE_FUNCTION real GetGamma(real P = 0.0, real rho = 0.0) const {return gamma;}
  void Refresh(DataBlock &, real) {}  // Refresh the eos (recompute coefficients and tables)

  KOKKOS_INLINE_FUNCTION
  real GetWaveSpeed(int k, int j, int i) const {
    Kokkos::abort("GetWaveSpeed should be used only for isothermal EOS");
    return 0;
  }
  KOKKOS_INLINE_FUNCTION
  real GetInternalEnergy(real P, real rho) const {
    return P/(gamma-1.0);
  }
  KOKKOS_INLINE_FUNCTION
  real GetPressure(real Eint, real rho) const {
    return Eint * (gamma-1.0);
  }

 private:
  real gamma;
};
#endif // FLUID_EOS_EOS_ADIABATIC_HPP_
