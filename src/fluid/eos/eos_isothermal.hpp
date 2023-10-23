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


// This is the isothermal implementation of the equation of state
class EquationOfState {
 public:
  EquationOfState() = default;

  EquationOfState(Input & input, DataBlock *data, std::string prefix) {
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

  void ShowConfig() {
    if(haveIsoSoundSpeed == Constant) {
    idfx::cout << "EquationOfState: isothermal with cs=" << isoSoundSpeed << "."
                << std::endl;
    } else if(haveIsoSoundSpeed == UserDefFunction) {
      idfx::cout << "EquationOfState: isothermal with user-defined cs function."
                  << std::endl;
      if(!isoSoundSpeedFunc) {
        IDEFIX_ERROR("No user-defined isothermal sound speed function has been enrolled.");
      }
    }
  }

  // In the isothermal EOS, gamma does not depend on the gas state,
  // So we add default values to 0 here so that GetGamma can be called without any argument
  KOKKOS_INLINE_FUNCTION real GetGamma(real P = 0.0, real rho = 0.0) const {
    return 1.0;
  }

  void Refresh(DataBlock &data, real t) {     // Refresh the coefficients (and tables)
  idfx::pushRegion("EquationOfState::Refresh");
    if(haveIsoSoundSpeed == UserDefFunction) {
      if(isoSoundSpeedFunc) {
        idfx::pushRegion("EquationOfState::UserDefSoundSpeed");
        isoSoundSpeedFunc(data, t, isoSoundSpeedArray);
        idfx::popRegion();
      } else {
        IDEFIX_ERROR("No user-defined isothermal sound speed function has been enrolled");
      }
    }
    idfx::popRegion();
  }

  KOKKOS_INLINE_FUNCTION real GetWaveSpeed(const int k, const int j, const int i) const {
    if(haveIsoSoundSpeed == UserDefFunction) {
      return isoSoundSpeedArray(k,j,i);
    } else {
      return isoSoundSpeed;
    }
  }
  KOKKOS_INLINE_FUNCTION real GetInternalEnergy(real P, real rho) const {
    Kokkos::abort("Internal energy is not defined in the isothermal EOS");
  }
  KOKKOS_INLINE_FUNCTION real GetPressure(real Eint, real rho) const {
    Kokkos::abort("Pressure is not defined in the isothermal EOS");
  }

  // Enroll user-defined isothermal sound speed
  void EnrollIsoSoundSpeed(IsoSoundSpeedFunc func) {
    if(this->haveIsoSoundSpeed != UserDefFunction) {
      IDEFIX_WARNING("Isothermal sound speed enrollment requires Hydro/csiso "
                  " to be set to userdef in .ini file");
    }
    this->isoSoundSpeedFunc = func;
  }

 private:
    real isoSoundSpeed;
    HydroModuleStatus haveIsoSoundSpeed{Disabled};
    IdefixArray3D<real> isoSoundSpeedArray;
    IsoSoundSpeedFunc isoSoundSpeedFunc{NULL};
};

#endif // FLUID_EOS_EOS_ISOTHERMAL_HPP_
