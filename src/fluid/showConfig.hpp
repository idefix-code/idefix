// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_SHOWCONFIG_HPP_
#define FLUID_SHOWCONFIG_HPP_

#include <string>

#include "idefix.hpp"
#include "fluid.hpp"
#include "dataBlock.hpp"

template<typename Phys>
void Fluid<Phys>::ShowConfig() {
  idfx::cout << Phys::prefix << ": ";
  if constexpr(!Phys::dust) {
    if constexpr(Phys::mhd) {
      idfx::cout << "solving MHD equations." << std::endl;
      #ifdef EVOLVE_VECTOR_POTENTIAL
        idfx::cout << Phys::prefix << ": Using EXPERIMENTAL vector potential formulation for MHD."
                  << std::endl;
      #endif
    } else {
      idfx::cout << "solving HD equations." << std::endl;
    }
  } else {
    idfx::cout << "solving pressure-less dust equations." << std::endl;
  }
  idfx::cout << Phys::prefix << ": Reconstruction: ";
  #if ORDER == 1
    idfx::cout << "1st order (donor cell)" << std::endl;
  #elif ORDER == 2
    idfx::cout << "2nd order (PLM Van Leer)" << std::endl;
  #elif ORDER == 3
    idfx::cout << "3rd order (LimO3)" << std::endl;
  #elif ORDER == 4
    idfx::cout << "4th order (PPM)" << std::endl;
  #endif



  if(haveRotation) {
    idfx::cout << Phys::prefix << ": Rotation ENABLED with Omega=" << this->OmegaZ << std::endl;
  }
  if(haveShearingBox) {
    idfx::cout << Phys::prefix << ": ShearingBox ENABLED with shear rate= " << this->sbS
               << " and Lx= " << sbLx << std::endl;
  }

  if(resistivityStatus.status != Disabled) {
    if(resistivityStatus.status == Constant) {
      idfx::cout << Phys::prefix << ": Ohmic resistivity ENABLED with constant resistivity eta="
                 << etaO << std::endl;
    } else if(resistivityStatus.status == UserDefFunction) {
      idfx::cout << Phys::prefix
                 << ": Ohmic resistivity ENABLED with user-defined resistivity function."
                 << std::endl;
      if(!ohmicDiffusivityFunc) {
        IDEFIX_ERROR("No user-defined Ihmic resistivity function has been enrolled.");
      }
    } else {
      IDEFIX_ERROR("Unknown Ohmic resistivity mode");
    }
    if(resistivityStatus.isExplicit) {
      idfx::cout << Phys::prefix << ": Ohmic resistivity uses an explicit time integration."
                 << std::endl;
    } else if(resistivityStatus.isRKL) {
      idfx::cout << Phys::prefix
                 << ": Ohmic resistivity uses a Runge-Kutta-Legendre time integration."
                 << std::endl;
    } else {
      IDEFIX_ERROR("Unknown time integrator for Ohmic resistivity");
    }
  }

  if(ambipolarStatus.status != Disabled) {
    if(ambipolarStatus.status == Constant) {
      idfx::cout << Phys::prefix << ": Ambipolar diffusion ENABLED with constant diffusivity xA="
                 << xA << std::endl;
    } else if(ambipolarStatus.status == UserDefFunction) {
      idfx::cout << Phys::prefix
                 << ": Ambipolar diffusion ENABLED with user-defined diffusivity function."
                 << std::endl;
      if(!ambipolarDiffusivityFunc) {
        IDEFIX_ERROR("No user-defined ambipolar diffusion function has been enrolled.");
      }
    } else {
      IDEFIX_ERROR("Unknown Ambipolar diffusion mode");
    }
    if(ambipolarStatus.isExplicit) {
      idfx::cout << Phys::prefix << ": Ambipolar diffusion uses an explicit time integration."
                 << std::endl;
    } else if(ambipolarStatus.isRKL) {
      idfx::cout << Phys::prefix
                 << ": Ambipolar diffusion uses a Runge-Kutta-Legendre time integration."
                 << std::endl;
    } else {
      IDEFIX_ERROR("Unknown time integrator for ambipolar diffusion");
    }
  }

  if(hallStatus.status != Disabled) {
    if(hallStatus.status == Constant) {
      idfx::cout << Phys::prefix << ": Hall effect ENABLED with constant diffusivity xH="
                 << xH << std::endl;
    } else if(hallStatus.status == UserDefFunction) {
      idfx::cout << Phys::prefix << ": Hall effect ENABLED with user-defined diffusivity function."
                 << std::endl;
      if(!hallDiffusivityFunc) {
        IDEFIX_ERROR("No user-defined Hall diffusivity function has been enrolled.");
      }
    } else {
      IDEFIX_ERROR("Unknown Hall effect mode");
    }
    if(hallStatus.isExplicit) {
      idfx::cout << Phys::prefix << ": Hall effect uses an explicit time integration." << std::endl;
    }  else {
      IDEFIX_ERROR("Unknown time integrator for Hall effect");
    }
  }

  if(haveTracer) {
    idfx::cout << Phys::prefix << ": " << this->nTracer << " tracers ENABLED for this fluid."
               << std::endl;
  }

  if(emfBoundaryFunc) {
    idfx::cout << Phys::prefix << ": user-defined EMF boundaries ENABLED." << std::endl;
  }
  if(userSourceTerm) {
    idfx::cout << Phys::prefix << ": user-defined source terms ENABLED." << std::endl;
  }

  if constexpr(Phys::eos) {
    eos->ShowConfig();
  }
  rSolver->ShowConfig();

  if constexpr(Phys::mhd) {
    emf->ShowConfig();
  }
  if(haveRKLParabolicTerms) {
    rkl->ShowConfig();
  }
  if(viscosityStatus.isExplicit || viscosityStatus.isRKL) {
    viscosity->ShowConfig();
  }
  if(thermalDiffusionStatus.status != Disabled) {
    thermalDiffusion->ShowConfig();
  }
  if(bragViscosityStatus.isExplicit || bragViscosityStatus.isRKL) {
    bragViscosity->ShowConfig();
  }
  if(bragThermalDiffusionStatus.status != Disabled) {
    bragThermalDiffusion->ShowConfig();
  }
  if(haveAxis) {
    boundary->axis->ShowConfig();
  }
  if(haveDrag) {
    drag->ShowConfig();
  }
}
#endif //FLUID_SHOWCONFIG_HPP_
