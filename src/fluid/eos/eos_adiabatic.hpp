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
  EquationOfState(Input & input, DataBlock *, std::string prefix) {
    this->gamma = input.GetOrSet<real>(prefix,"gamma",0, 5.0/3.0);
    this->prefix = prefix;
  }
    
  void ShowConfig() {
    idfx::cout << prefix << ": EOS: ideal with gamma=" << this->gamma << std::endl;
  };

  KOKKOS_INLINE_FUNCTION real GetGamma(real P, real rho) const {return gamma;}
  void Refresh(DataBlock &, real) {}  // Refresh the eos (recompute coefficients and tables)

  KOKKOS_INLINE_FUNCTION 
  real GetWaveSpeed(real P, real rho) const {
    return std::sqrt(gamma*(P/rho));
  }
  KOKKOS_INLINE_FUNCTION 
  real GetInternalEnergy(real P, real rho) const {
    return P/(gamma-1.0);
  };
  KOKKOS_INLINE_FUNCTION 
  real GetPressure(real Eint, real rho) const {
    return Eint * (gamma-1.0);
  };

 private:
  real gamma;
  std::string prefix;

};
#endif // FLUID_EOS_ADIABATIC_HPP