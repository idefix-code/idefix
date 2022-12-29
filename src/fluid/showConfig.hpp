// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>

#include "idefix.hpp"
#include "fluid.hpp"
#include "dataBlock.hpp"

template<typename Phys>
void Fluid<Phys>::ShowConfig() {
  idfx::cout << "Fluid: ";
  #if MHD == YES
    idfx::cout << "solving MHD equations." << std::endl;
    #ifdef EVOLVE_VECTOR_POTENTIAL
      idfx::cout << "Fluid: Using EXPERIMENTAL vector potential formulation for MHD." << std::endl;
    #endif
  #else
    idfx::cout << "solving HD equations." << std::endl;
  #endif
  idfx::cout << "Fluid: Reconstruction: ";
  #if ORDER == 1
    idfx::cout << "1st order (donor cell)" << std::endl;
  #elif ORDER == 2
    idfx::cout << "2nd order (PLM Van Leer)" << std::endl;
  #elif ORDER == 3
    idfx::cout << "3rd order (LimO3)" << std::endl;
  #elif ORDER == 4
    idfx::cout << "4th order (PPM)" << std::endl;
  #endif

  #if HAVE_ENERGY
    idfx::cout << "Fluid: EOS: ideal with gamma=" << this->gamma << std::endl;
  #endif
  #ifdef ISOTHERMAL
    if(haveIsoSoundSpeed == Constant) {
      idfx::cout << "Fluid: EOS: isothermal with cs=" << isoSoundSpeed << "." << std::endl;
    } else if(haveIsoSoundSpeed == UserDefFunction) {
      idfx::cout << "Fluid: EOS: isothermal with user-defined cs function." << std::endl;
      if(!isoSoundSpeedFunc) {
        IDEFIX_ERROR("No user-defined isothermal sound speed function has been enrolled.");
      }
    }
  #endif// ISOTHERMAL
  idfx::cout << "Fluid: Riemann solver: ";
  switch(mySolver) {
    case TVDLF:
      idfx::cout << "tvdlf." << std::endl;
      break;
    case HLL:
      idfx::cout << "hll." << std::endl;
      break;
    #if MHD==YES
      case HLLD:
        idfx::cout << "hlld." << std::endl;
        break;
    #else
      case HLLC:
        idfx::cout << "hllc." << std::endl;
        break;
    #endif
    case ROE:
      idfx::cout << "roe." << std::endl;
      break;
    default:
      IDEFIX_ERROR("Unknown Riemann solver");
  }
  if(haveRotation) {
    idfx::cout << "Fluid: Rotation ENABLED with Omega=" << this->OmegaZ << std::endl;
  }
  if(haveShearingBox) {
    idfx::cout << "Fluid: ShearingBox ENABLED with shear rate= " << this->sbS
               << " and Lx= " << sbLx << std::endl;
  }
  if(haveShockFlattening) {
    idfx::cout << "Fluid: Shock Flattening ENABLED." << std::endl;
  }


  if(resistivityStatus.status != Disabled) {
    if(resistivityStatus.status == Constant) {
      idfx::cout << "Fluid: Ohmic resistivity ENABLED with constant resistivity eta="
                 << etaO << std::endl;
    } else if(resistivityStatus.status == UserDefFunction) {
      idfx::cout << "Fluid: Ohmic resistivity ENABLED with user-defined resistivity function."
                 << std::endl;
      if(!ohmicDiffusivityFunc) {
        IDEFIX_ERROR("No user-defined Ihmic resistivity function has been enrolled.");
      }
    } else {
      IDEFIX_ERROR("Unknown Ohmic resistivity mode");
    }
    if(resistivityStatus.isExplicit) {
      idfx::cout << "Fluid: Ohmic resistivity uses an explicit time integration." << std::endl;
    } else if(resistivityStatus.isRKL) {
      idfx::cout << "Fluid: Ohmic resistivity uses a Runge-Kutta-Legendre time integration."
                 << std::endl;
    } else {
      IDEFIX_ERROR("Unknown time integrator for Ohmic resistivity");
    }
  }

  if(ambipolarStatus.status != Disabled) {
    if(ambipolarStatus.status == Constant) {
      idfx::cout << "Fluid: Ambipolar diffusion ENABLED with constant diffusivity xA="
                 << xA << std::endl;
    } else if(ambipolarStatus.status == UserDefFunction) {
      idfx::cout << "Fluid: Ambipolar diffusion ENABLED with user-defined diffusivity function."
                 << std::endl;
      if(!ambipolarDiffusivityFunc) {
        IDEFIX_ERROR("No user-defined ambipolar diffusion function has been enrolled.");
      }
    } else {
      IDEFIX_ERROR("Unknown Ambipolar diffusion mode");
    }
    if(ambipolarStatus.isExplicit) {
      idfx::cout << "Fluid: Ambipolar diffusion uses an explicit time integration." << std::endl;
    } else if(ambipolarStatus.isRKL) {
      idfx::cout << "Fluid: Ambipolar diffusion uses a Runge-Kutta-Legendre time integration."
                 << std::endl;
    } else {
      IDEFIX_ERROR("Unknown time integrator for ambipolar diffusion");
    }
  }

  if(hallStatus.status != Disabled) {
    if(hallStatus.status == Constant) {
      idfx::cout << "Fluid: Hall effect ENABLED with constant diffusivity xH="
                 << xH << std::endl;
    } else if(hallStatus.status == UserDefFunction) {
      idfx::cout << "Fluid: Hall effect ENABLED with user-defined diffusivity function."
                 << std::endl;
      if(!hallDiffusivityFunc) {
        IDEFIX_ERROR("No user-defined Hall diffusivity function has been enrolled.");
      }
    } else {
      IDEFIX_ERROR("Unknown Hall effect mode");
    }
    if(hallStatus.isExplicit) {
      idfx::cout << "Fluid: Hall effect uses an explicit time integration." << std::endl;
    }  else {
      IDEFIX_ERROR("Unknown time integrator for Hall effect");
    }
  }

  if(emfBoundaryFunc) {
    idfx::cout << "Fluid: user-defined EMF boundaries ENABLED." << std::endl;
  }
  if(userSourceTerm) {
    idfx::cout << "Fluid: user-defined source terms ENABLED." << std::endl;
  }
  #if MHD == YES
    emf->ShowConfig();
  #endif
  if(haveRKLParabolicTerms) {
    rkl->ShowConfig();
  }
  if(viscosityStatus.isExplicit || viscosityStatus.isRKL) {
    viscosity.ShowConfig();
  }
  if(thermalDiffusionStatus.isExplicit || thermalDiffusionStatus.isRKL) {
    thermalDiffusion.ShowConfig();
  }
  if(haveAxis) {
    myAxis->ShowConfig();
  }
}