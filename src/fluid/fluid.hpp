// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_FLUID_HPP_
#define FLUID_FLUID_HPP_

#include <string>
#include <vector>
#include "idefix.hpp"

#include "grid.hpp"
#include "fluid_defs.hpp"
#include "electroMotiveForce.hpp"
#include "viscosity.hpp"
#include "thermalDiffusion.hpp"
#include "axis.hpp"
#include "hydroboundary.hpp"
#include "shockFlattening.hpp"
#include "selfGravity.hpp"

// forward class declaration
class DataBlock;

template<typename Phys>
class Fluid {
 public:
  Fluid();
  void Init(Input &, Grid &, DataBlock *);
  void ConvertConsToPrim();
  void ConvertPrimToCons();
  template <int> void CalcRiemannFlux(const real);
  template <int> void CalcParabolicFlux(const real);
  template <int> void AddNonIdealMHDFlux(const real);
  template <int> void CalcRightHandSide(real, real );
  void CalcCurrent();
  void AddSourceTerms(real, real );
  void CoarsenFlow(IdefixArray4D<real>&);
  void CoarsenMagField(IdefixArray4D<real>&);
  void CoarsenVectorPotential();
  real GetGamma();
  real CheckDivB();
  void ResetStage();
  void ShowConfig();

  // Our boundary conditions
  HydroBoundary boundary;

  // Source terms
  bool haveSourceTerms{false};

  // Parabolic terms
  bool haveExplicitParabolicTerms{false};
  bool haveRKLParabolicTerms{false};

  // Current
  bool haveCurrent{false};
  bool needExplicitCurrent{false};
  bool needRKLCurrent{false};

  // Nonideal MHD effects coefficients
  ParabolicModuleStatus resistivityStatus, ambipolarStatus, hallStatus;

  // Whether or not we have viscosity
  ParabolicModuleStatus viscosityStatus;

  // Whether or not we have thermal diffusion
  ParabolicModuleStatus thermalDiffusionStatus;

  // Viscosity object
  Viscosity viscosity;

  // Thermal Diffusion object
  ThermalDiffusion thermalDiffusion;

  // Whether or not we have to treat the axis
  bool haveAxis{false};
  Axis myAxis;

  // Rotation vector
  bool haveRotation{false};
  real OmegaZ;

  bool haveShearingBox{false};
  // Shear rate for shearing box problems
  real sbS;
  // Box width for shearing box problems
  real sbLx;

  // ShockFlattening
  bool haveShockFlattening{false};
  ShockFlattening shockFlattening;

  // Enroll user-defined boundary conditions
  void EnrollUserDefBoundary(UserDefBoundaryFunc);
  void EnrollInternalBoundary(InternalBoundaryFunc);
  void EnrollEmfBoundary(EmfBoundaryFunc);
  void EnrollFluxBoundary(UserDefBoundaryFunc);

  // Add some user source terms
  void EnrollUserSourceTerm(SrcTermFunc);

  // DEPRECATED gravity enrollment
  void EnrollGravPotential(GravPotentialFunc);
  void EnrollBodyForce(BodyForceFunc);

  // Enroll user-defined ohmic, ambipolar and Hall diffusivities
  void EnrollOhmicDiffusivity(DiffusivityFunc);
  void EnrollAmbipolarDiffusivity(DiffusivityFunc);
  void EnrollHallDiffusivity(DiffusivityFunc);

  // Enroll user-defined isothermal sound speed
  void EnrollIsoSoundSpeed(IsoSoundSpeedFunc);

  // Riemann Solvers
#if MHD == YES
  template<const int>
    void HlldMHD();
  template<const int>
    void HllMHD();
  template<const int>
    void RoeMHD();
  template<const int>
    void TvdlfMHD();
#else
  template<const int>
    void HllcHD();
  template<const int>
    void HllHD();
  template<const int>
    void RoeHD();
  template<const int>
    void TvdlfHD();
#endif

  // Arrays required by the Hydro object
  IdefixArray4D<real> Vc;      // Main cell-centered primitive variables index
  IdefixArray4D<real> Vs;      // Main face-centered varariables
  IdefixArray4D<real> Ve;      // Main edge-centered varariables (only when EVOLVE_VECTOR_POTENTIAL)
  IdefixArray4D<real> Uc;      // Main cell-centered conservative variables
  IdefixArray4D<real> J;       // Electrical current
                               // (only defined when non-ideal MHD effects are enabled)

  // Name of the fields (used in outputs)
  std::vector<std::string> VcName;
  std::vector<std::string> VsName;
  std::vector<std::string> VeName;

  // Storing all of the electromotive forces
  ElectroMotiveForce emf;

  // Required by time integrator
  IdefixArray4D<real> Uc0;
  IdefixArray4D<real> Vs0;
  IdefixArray4D<real> Ve0;
  IdefixArray3D<real> InvDt;

  IdefixArray4D<real> FluxRiemann;
  IdefixArray3D<real> dMax;    // Maximum diffusion speed



 private:
  friend class ElectroMotiveForce;
  friend class Viscosity;
  friend class ThermalDiffusion;
  friend class Fargo;
  friend class Axis;
  friend class RKLegendre;
  friend class HydroBoundary;
  friend class ShockFlattening;

  // Isothermal EOS parameters
  real isoSoundSpeed;
  HydroModuleStatus haveIsoSoundSpeed{Disabled};
  IdefixArray3D<real> isoSoundSpeedArray;
  IsoSoundSpeedFunc isoSoundSpeedFunc{NULL};

  // Adiabatic EOS parameters
  real gamma;

  Solver mySolver;

  DataBlock *data;

  // Emf boundary conditions
  bool haveEmfBoundary{false};
  EmfBoundaryFunc emfBoundaryFunc{NULL};

  // User defined source term
  SrcTermFunc userSourceTerm{NULL};
  bool haveUserSourceTerm{false};

  real etaO, xH, xA;  // Ohmic resistivity, Hall, ambipolar (when constant)

  // Ohmic, Hall and ambipolar diffusivity (when function-defined)
  DiffusivityFunc ohmicDiffusivityFunc{NULL};
  DiffusivityFunc ambipolarDiffusivityFunc{NULL};
  DiffusivityFunc hallDiffusivityFunc{NULL};

  IdefixArray3D<real> cMax;    // Maximum propagation speed

  // Nonideal effect diffusion coefficient (only allocated when needed)
  IdefixArray3D<real> etaOhmic;
  IdefixArray3D<real> xHall;
  IdefixArray3D<real> xAmbipolar;
};

#include "../physics.hpp"

using Hydro = Fluid<Physics>;

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
    emf.ShowConfig();
  #endif
  if(viscosityStatus.isExplicit || viscosityStatus.isRKL) {
    viscosity.ShowConfig();
  }
  if(thermalDiffusionStatus.isExplicit || thermalDiffusionStatus.isRKL) {
    thermalDiffusion.ShowConfig();
  }
  if(haveAxis) {
    myAxis.ShowConfig();
  }
}

#include "addSourceTerms.hpp"
#include "enroll.hpp"
#include "init.hpp"
#include "calcCurrent.hpp"
#include "coarsenFlow.hpp"
#include "convertConsToPrim.hpp"
#endif // FLUID_FLUID_HPP_
