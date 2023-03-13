// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_ENROLL_HPP_
#define FLUID_ENROLL_HPP_

#include "dataBlock.hpp"
template<typename Phys>
void Fluid<Phys>::EnrollIsoSoundSpeed(IsoSoundSpeedFunc myFunc) {
  if(this->haveIsoSoundSpeed != UserDefFunction) {
    IDEFIX_WARNING("Isothermal sound speed enrollment requires Hydro/csiso "
                 " to be set to userdef in .ini file");
  }
  if constexpr(!Phys::isothermal) {
    IDEFIX_ERROR("Isothermal sound speed enrollment requires ISOTHERMAL to be defined in"
                 "definitions.hpp");
  }
  this->isoSoundSpeedFunc = myFunc;
}

template<typename Phys>
void Fluid<Phys>::EnrollUserDefBoundary(UserDefBoundaryFunc myFunc) {
  // This is a proxy for userdef enrollment
  boundary->EnrollUserDefBoundary(myFunc);
}

template<typename Phys>
void Fluid<Phys>::EnrollFluxBoundary(UserDefBoundaryFunc myFunc) {
  // This is a proxy for userdef enrollment
  boundary->EnrollFluxBoundary(myFunc);
}

template<typename Phys>
void Fluid<Phys>::EnrollInternalBoundary(InternalBoundaryFunc myFunc) {
  // This is a proxy for userdef enrollment
  boundary->EnrollInternalBoundary(myFunc);
}

template<typename Phys>
void Fluid<Phys>::EnrollUserSourceTerm(SrcTermFunc myFunc) {
  this->userSourceTerm = myFunc;
  this->haveUserSourceTerm = true;
  this->haveSourceTerms = true;
}


template<typename Phys>
void Fluid<Phys>::EnrollEmfBoundary(EmfBoundaryFunc myFunc) {
  if constexpr(!Phys::mhd) {
    IDEFIX_ERROR("This function can only be used with the MHD solver.");
  }
  this->emfBoundaryFunc = myFunc;
  this->haveEmfBoundary = true;
}

template<typename Phys>
void Fluid<Phys>::EnrollOhmicDiffusivity(DiffusivityFunc myFunc) {
  if constexpr(!Phys::mhd) {
    IDEFIX_ERROR("This function can only be used with the MHD solver.");
  }
  if(this->resistivityStatus.status < UserDefFunction) {
    IDEFIX_WARNING("Ohmic diffusivity enrollment requires Hydro/Resistivity "
                 "to be set to userdef in .ini file");
  }
  this->ohmicDiffusivityFunc = myFunc;
}

template<typename Phys>
void Fluid<Phys>::EnrollAmbipolarDiffusivity(DiffusivityFunc myFunc) {
  if constexpr(!Phys::mhd) {
    IDEFIX_ERROR("This function can only be used with the MHD solver.");
  }
  if(this->ambipolarStatus.status < UserDefFunction) {
    IDEFIX_WARNING("Ambipolar diffusivity enrollment requires Hydro/Ambipolar "
                 "to be set to userdef in .ini file");
  }
  this->ambipolarDiffusivityFunc = myFunc;
}

template<typename Phys>
void Fluid<Phys>::EnrollHallDiffusivity(DiffusivityFunc myFunc) {
  if constexpr(!Phys::mhd) {
    IDEFIX_ERROR("This function can only be used with the MHD solver.");
  }
  if(this->hallStatus.status < UserDefFunction) {
    IDEFIX_WARNING("Hall diffusivity enrollment requires Hydro/Hall "
                 "to be set to userdef in .ini file");
  }
  this->hallDiffusivityFunc = myFunc;
}

template<typename Phys>
real Fluid<Phys>::GetGamma() {
  return(this->gamma);
}

template<typename Phys>
void Fluid<Phys>::ResetStage() {
  // Reset variables required at the beginning of each stage
  // (essentially linked to timestep evaluation)
  idfx::pushRegion("Fluid::ResetStage");

  IdefixArray3D<real> InvDt=this->InvDt;

  idefix_for("HydroResetStage",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      InvDt(k,j,i) = ZERO_F;
  });

  idfx::popRegion();
}

#endif //FLUID_ENROLL_HPP_
