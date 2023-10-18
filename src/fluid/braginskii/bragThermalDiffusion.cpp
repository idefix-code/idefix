// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// This source code is largely inspired from the th_flux of Pluto4.2
// ((c) P. Tzeferacos & A. Mignone)

// Implementation of monotonicity-preserving heat flux following Sharma & Hammett 2007,
// Jour. of Comp. Physics

#include <string>

#include "bragThermalDiffusion.hpp"
#include "dataBlock.hpp"
#include "fluid.hpp"
#include "eos.hpp"



void BragThermalDiffusion::ShowConfig() {
  if(status.status==Constant) {
    idfx::cout << "Braginskii Thermal Diffusion: ENABLED with constant diffusivity kpar="
                    << this->kpar << " and knor=" << this->knor << " ."<< std::endl;
  } else if (status.status==UserDefFunction) {
    idfx::cout << "Braginskii Thermal Diffusion: ENABLED with user-defined diffusivity function."
                   << std::endl;
    if(!diffusivityFunc) {
      IDEFIX_ERROR("No braginskii thermal diffusion function has been enrolled");
    }
  } else {
    IDEFIX_ERROR("Unknown braginskii thermal diffusion mode");
  }
  if(status.isExplicit) {
    idfx::cout << "Braginskii Thermal Diffusion: uses an explicit time integration." << std::endl;
  } else if(status.isRKL) {
    idfx::cout << "Braginskii Thermal Diffusion: uses a Runge-Kutta-Legendre time integration."
                << std::endl;
  } else {
    IDEFIX_ERROR("Unknown time integrator for braginskii thermal diffusion.");
  }
  if(haveSlopeLimiter) {
    idfx::cout << "Braginskii Thermal Diffusion: uses a slope limiter." << std::endl;
  }
}

void BragThermalDiffusion::EnrollBragThermalDiffusivity(BragDiffusivityFunc myFunc) {
  if(this->status.status != UserDefFunction) {
    IDEFIX_WARNING("Braginskii thermal diffusivity enrollment requires Hydro/BragThermalDiffusion "
                 "to be set to userdef in .ini file");
  }
  this->diffusivityFunc = myFunc;
}

void BragThermalDiffusion::AddBragDiffusiveFlux(int dir, const real t,
                                                const IdefixArray4D<real> &Flux) {
  idfx::pushRegion("BragThermalDiffusion::AddBragDiffusiveFlux");
  switch(limiter) {
    case PLMLimiter::VanLeer:
      this->AddBragDiffusiveFluxLim<PLMLimiter::VanLeer>(dir,t,Flux);
      break;
    case PLMLimiter::McLim:
      this->AddBragDiffusiveFluxLim<PLMLimiter::McLim>(dir,t,Flux);
      break;
    case PLMLimiter::MinMod:
      this->AddBragDiffusiveFluxLim<PLMLimiter::MinMod>(dir,t,Flux);
      break;
    default:
      IDEFIX_ERROR("The slope limiter for the Braginskii heat flux is not defined.");
      break;
  }
  idfx::popRegion();
}
