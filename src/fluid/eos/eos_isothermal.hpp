// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_EOS_EOS_ISOTHERMAL_HPP_
#define FLUID_EOS_EOS_ISOTHERMAL_HPP_

#include <string>
#include "idefix.hpp"
#include "input.hpp"
#include "fluid_defs.hpp"
#include "dataBlock.hpp"


// This is the adiabatic (constant gamma) implementation of the equation of state
class EquationOfState {
 public:
  EquationOfState(Input &, DataBlock *, std::string);

  void ShowConfig();

  real GetGamma(real P, real rho);            // The gamma of this EOS
  void Refresh(DataBlock &);  // Refresh the eos (recompute coefficients and tables)
  KOKKOS_INLINE_FUNCTION real GetWaveSpeed(real P, real rho);
  KOKKOS_INLINE_FUNCTION real GetWaveSpeed(int i, int j, int k) {
    if(haveIsoSoundSpeed == UserDefFunction) {
      return isoSoundSpeedArray(k,j,i);
    } else {
      return isoSoundSpeed;
    }
  }
  KOKKOS_INLINE_FUNCTION real GetInternalEnergy(real P, real rho);
  KOKKOS_INLINE_FUNCTION real GetPressure(real Eint, real rho);

  // Enroll user-defined isothermal sound speed
  void EnrollIsoSoundSpeed(IsoSoundSpeedFunc);

 private:
    real isoSoundSpeed;
    HydroModuleStatus haveIsoSoundSpeed{Disabled};
    IdefixArray3D<real> isoSoundSpeedArray;
    IsoSoundSpeedFunc isoSoundSpeedFunc{NULL};
    std::string prefix;

};

EquationOfState::EquationOfState(Input & input, DataBlock *data, std::string prefix) {
  this->prefix = prefix;
  std::string isoString = input.Get<std::string>(prefix,"csiso",0);
  if(isoString.compare("constant") == 0) {
    this->haveIsoSoundSpeed = Constant;
    this->isoSoundSpeed = input.Get<real>(prefix,"csiso",1);
  } else if(isoString.compare("userdef") == 0) {
    this->haveIsoSoundSpeed = UserDefFunction;
    this->isoSoundSpeedArray = IdefixArray3D<real>(prefix+"_csIso",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  } else {
    IDEFIX_ERROR("csiso admits only constant or userdef entries");
  }
}

void EquationOfState::Refresh(DataBlock &data, real t) {
  idfx::pushRegion("EquationOfState::Refresh");
  if(haveIsoSoundSpeed == UserDefFunction) {
    if(isoSoundSpeedFunc)
      idfx::pushRegion("EquationOfState::UserDefSoundSpeed");
      isoSoundSpeedFunc(data, t, isoSoundSpeedArray);
      idfx::popRegion();
    else
      IDEFIX_ERROR("No user-defined isothermal sound speed function has been enrolled");
  }
  idfx::popRegion();
}

void EquationOfState::EnrollIsoSoundSpeed(IsoSoundSpeedFunc func) {
  if(this->haveIsoSoundSpeed != UserDefFunction) {
    IDEFIX_WARNING("Isothermal sound speed enrollment requires Hydro/csiso "
                 " to be set to userdef in .ini file");
  }
  this->isoSoundSpeedFunc = func;
}

void EquationOfState::ShowConfig() {
  if(haveIsoSoundSpeed == Constant) {
    idfx::cout << prefix << ": EOS: isothermal with cs=" << isoSoundSpeed << "."
                << std::endl;
  } else if(haveIsoSoundSpeed == UserDefFunction) {
    idfx::cout << prefix << ": EOS: isothermal with user-defined cs function."
                << std::endl;
    if(!isoSoundSpeedFunc) {
      IDEFIX_ERROR("No user-defined isothermal sound speed function has been enrolled.");
    }
  }
}

#endif // FLUID_EOS_ISOTHERMAL_HPP