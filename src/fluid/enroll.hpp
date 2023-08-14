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
  if constexpr(!Phys::isothermal) {
    IDEFIX_ERROR("Isothermal sound speed enrollment requires ISOTHERMAL to be defined in"
                 "definitions.hpp");
  } else {
    #ifdef ISOTHERMAL
    eos->EnrollIsoSoundSpeed(myFunc);
    #endif
  }
}

template<typename Phys>
template<typename T>
void Fluid<Phys>::EnrollUserDefBoundary(T myFunc) {
  // This is a proxy for userdef enrollment
  boundary->EnrollUserDefBoundary(myFunc);
}

template<typename Phys>
template<typename T>
void Fluid<Phys>::EnrollInternalBoundary(T myFunc) {
  // This is a proxy for internal boundary enrollment
  boundary->EnrollInternalBoundary(myFunc);
}

template<typename Phys>
template<typename T>
void Fluid<Phys>::EnrollFluxBoundary(T myFunc) {
  // This is a proxy for userdef enrollment
  boundary->EnrollFluxBoundary(myFunc);
}

template<typename Phys>
void Fluid<Phys>::EnrollUserSourceTerm(SrcTermFunc<Phys> myFunc) {
  this->userSourceTerm = myFunc;
  this->haveUserSourceTerm = true;
  this->haveSourceTerms = true;
}

// Deprecated enrollment function
template<typename Phys>
void Fluid<Phys>::EnrollUserSourceTerm(SrcTermFuncOld myFunc) {
  std::stringstream msg;
  msg << "The old signature for user-defined source terms " << std::endl
      << "(DataBlock &, const real t, const real dt)" << std::endl
      << "is deprecated. You should now use "<< std::endl
      << "(Fluid<Phys> *,  const real t, const real dt)" << std::endl
      << "With the Phys of your choice (DefaultPhysics, DustPhysics...)" << std::endl;

  IDEFIX_WARNING(msg);

  this->userSourceTermOld = myFunc;
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
