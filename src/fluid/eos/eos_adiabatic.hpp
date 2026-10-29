// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/fluid/eos/eos_adiabatic.hpp
//
// Last modified : 10/2023
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2023)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_EOS_EOS_ADIABATIC_HPP_
#define FLUID_EOS_EOS_ADIABATIC_HPP_

#include <string>
#include "idefix.hpp"
#include "input.hpp"

// This is the adiabatic (constant gamma) implementation
// of the equation of state (device-compatible subset)
class EquationOfStateDevice {
 public:
  EquationOfStateDevice() = default;

  EquationOfStateDevice(Input & input, DataBlock *, std::string prefix) {
    this->gamma = input.GetOrSet<real>(prefix,"gamma",0, 5.0/3.0);
  }

  // First adiabatic exponent. In the ideal EOS, gamma does not depend on the gas state,
  // So we add default values to 0 here so that GetGamma can be called without any argument
  KOKKOS_INLINE_FUNCTION real GetGamma(real P = 0.0, real rho = 0.0) const {return gamma;}

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

 protected:
  real gamma;
};

// This is the adiabatic (constant gamma) implementation of the equation of state (full class)
class EquationOfState : public EquationOfStateDevice {
 public:
  EquationOfState() = default;

  EquationOfState(Input & input, DataBlock *data, std::string prefix)
    : EquationOfStateDevice(input, data, prefix) {}

  void ShowConfig() {
    idfx::cout << "EquationOfState: ideal with gamma=" << this->gamma << std::endl;
  }

  void Refresh(DataBlock &, real) {}  // Refresh the eos (recompute coefficients and tables)
};
#endif // FLUID_EOS_EOS_ADIABATIC_HPP_
