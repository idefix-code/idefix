// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_EOS_EOS_TEMPLATE_HPP_
#define FLUID_EOS_EOS_TEMPLATE_HPP_

#include <string>
#include "idefix.hpp"
#include "input.hpp"

// This is a template for a generalized custom equation of state
class EquationOfState {
 public:
  EquationOfState() = default;

  // Add here the required steps to initialize your equation of state
  EquationOfState(Input & input, DataBlock *, std::string prefix) {
    // ....
  }

  // Information message displayed before entering the main loop
  void ShowConfig() {
    idfx::cout << "EquationOfState: Custom." << std::endl;
  }

  // First adiabatic exponent.
  KOKKOS_INLINE_FUNCTION real GetGamma(real P , real rho ) const {
    real gamma;
    // Compute gamma (needed for sound speed estimations, 5/3 would work in general)
    return gamma;
  }

  // Refresh the eos (recompute coefficients and tables)
  void Refresh(DataBlock &data, real t) {
    // ....
  }

  // This function is used only when the isothermal approximation is enabled. Not needed here
  KOKKOS_INLINE_FUNCTION
  real GetWaveSpeed(int k, int j, int i) const {
    Kokkos::abort("GetWaveSpeed should be used only for isothermal EOS");
    return 0;
  }

  // Compute the internal energy from pressure and density
  KOKKOS_INLINE_FUNCTION
  real GetInternalEnergy(real P, real rho) const {
    real eint; // = ...
    return eint;
  }

  // Compute the pressure from internal energy and density
  KOKKOS_INLINE_FUNCTION
  real GetPressure(real Eint, real rho) const {
    real P; // = ...
    return P;
  }

 private:
  real gamma;
  // Add here the internal variables required by your equation of state.
};


#endif // FLUID_EOS_EOS_TEMPLATE_HPP_
