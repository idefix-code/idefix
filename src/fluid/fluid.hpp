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
#include "output.hpp"
#include "fluid_defs.hpp"
#include "thermalDiffusion.hpp"
#include "selfGravity.hpp"
#include "vtk.hpp"
#include "dump.hpp"

// forward class declaration
class DataBlock;
template<typename Phys>
class Boundary;

template<typename Phys>
class Axis;

template<typename Phys>
class ConstrainedTransport;

template<typename Phys>
class RKLegendre;

template<typename Phys>
class RiemannSolver;

class Viscosity;

template<typename Phys>
class Fluid {
 public:
  Fluid( Grid &, Input&, DataBlock *);
  void ConvertConsToPrim();
  void ConvertPrimToCons();
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
  void EvolveStage(const real, const real);
  void ResetStage();
  void ShowConfig();

  // Our boundary conditions
  std::unique_ptr<Boundary<Phys>> boundary;

  // Source terms
  bool haveSourceTerms{false};

  // Parabolic terms
  bool haveExplicitParabolicTerms{false};
  bool haveRKLParabolicTerms{false};

  std::unique_ptr<RKLegendre<Phys>> rkl;

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
  std::unique_ptr<Viscosity> viscosity;

  // Thermal Diffusion object
  ThermalDiffusion thermalDiffusion;

  // Whether or not we have to treat the axis
  bool haveAxis{false};
  std::unique_ptr<Axis<Phys>> myAxis;

  // Rotation vector
  bool haveRotation{false};
  real OmegaZ;

  bool haveShearingBox{false};
  // Shear rate for shearing box problems
  real sbS;
  // Box width for shearing box problems
  real sbLx;


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
  std::unique_ptr<ConstrainedTransport<Phys>> emf;

  // Required by time integrator
  IdefixArray3D<real> InvDt;

  IdefixArray4D<real> FluxRiemann;
  IdefixArray3D<real> dMax;    // Maximum diffusion speed

  std::unique_ptr<RiemannSolver<Phys>> rSolver;


 private:
  friend class ConstrainedTransport<Phys>;
  friend class Viscosity;
  friend class ThermalDiffusion;
  friend class Fargo;
  friend class Axis<Phys>;
  friend class RKLegendre<Phys>;
  friend class Boundary<Phys>;
  friend class ShockFlattening;
  friend class RiemannSolver<Phys>;

  // Isothermal EOS parameters
  real isoSoundSpeed;
  HydroModuleStatus haveIsoSoundSpeed{Disabled};
  IdefixArray3D<real> isoSoundSpeedArray;
  IsoSoundSpeedFunc isoSoundSpeedFunc{NULL};

  // Adiabatic EOS parameters
  real gamma;

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

  // Loop on dimensions
  template <int dir>
  void LoopDir(const real, const real);
};

#include "../physics.hpp"
#include "dataBlock.hpp"
#include "boundary.hpp"
#include "constrainedTransport.hpp"
#include "axis.hpp"
#include "rkl.hpp"
#include "riemannSolver.hpp"
#include "viscosity.hpp"

using Hydro = Fluid<Physics>;


template<typename Phys>
Fluid<Phys>::Fluid(Grid &grid, Input &input, DataBlock *datain) {
  idfx::pushRegion("Fluid::Init");
  // Save the datablock to which we are attached from now on
  this->data = datain;

  #if ORDER < 1 || ORDER > 4
     IDEFIX_ERROR("Reconstruction at chosen order is not implemented. Check your definitions file");
  #endif


  #if HAVE_ENERGY
    this->gamma = input.GetOrSet<real>(std::string(Phys::prefix),"gamma",0, 5.0/3.0);
  #endif

  #ifdef ISOTHERMAL
    std::string isoString = input.Get<std::string>(std::string(Phys::prefix),"csiso",0);
    if(isoString.compare("constant") == 0) {
      this->haveIsoSoundSpeed = Constant;
      this->isoSoundSpeed = input.Get<real>(std::string(Phys::prefix),"csiso",1);
    } else if(isoString.compare("userdef") == 0) {
      this->haveIsoSoundSpeed = UserDefFunction;
    } else {
      IDEFIX_ERROR("csiso admits only constant or userdef entries");
    }
  #else
    // set the isothermal soundspeed, even though it will not be used
    this->isoSoundSpeed = -1.0;
  #endif

  // Source terms (always activated when non-cartesian geometry because of curvature source terms)
#if GEOMETRY == CARTESIAN
  this->haveSourceTerms = false;
#else
  this->haveSourceTerms = true;
#endif

  // Check whether we have rotation
  int rotation = input.CheckEntry(std::string(Phys::prefix),"rotation");

  if(rotation>=0 ) {
    this->haveSourceTerms = true;
    this->haveRotation = true;
    if(rotation >1 ) IDEFIX_ERROR("Rotation needs a a single rotation velocity (Omega_Z) \
                                   in idefix.ini");
    this->OmegaZ = input.Get<real>(std::string(Phys::prefix),"rotation",0);
  }

  // Check whether we have shearing box
  int shearingbox = input.CheckEntry(std::string(Phys::prefix),"shearingBox");

  if(shearingbox>=0 ) {
    this->haveShearingBox = true;
    this->haveSourceTerms = true;
    if(shearingbox != 1) {
      IDEFIX_ERROR("Shearing box needs a scalar value for the shear rate in idefix.ini");
    }

    this->sbS = input.Get<real>(std::string(Phys::prefix),"shearingBox",0);
    // Get box size
    this->sbLx = grid.xend[IDIR] - grid.xbeg[IDIR];
  }


  ///////////////////////
  // Parabolic terms
  ///////////////////////


  // Check whether viscosity is enabled, if so, init a viscosity object
  if(input.CheckEntry(std::string(Phys::prefix),"viscosity")>=0) {
    std::string opType = input.Get<std::string>(std::string(Phys::prefix),"viscosity",0);
    if(opType.compare("explicit") == 0 ) {
      haveExplicitParabolicTerms = true;
      viscosityStatus.isExplicit = true;
    } else if(opType.compare("rkl") == 0 ) {
      haveRKLParabolicTerms = true;
      viscosityStatus.isRKL = true;
    } else {
      std::stringstream msg;
      msg  << "Unknown integration type for viscosity: " << opType;
      IDEFIX_ERROR(msg);
    }
    this->viscosity = std::unique_ptr<Viscosity>(new Viscosity(input, grid, this));
  }

  // Check whether thermal diffusion is enabled, if so, init a thermal diffusion object
  if(input.CheckEntry(std::string(Phys::prefix),"TDiffusion")>=0) {
    std::string opType = input.Get<std::string>(std::string(Phys::prefix),"TDiffusion",0);
    if(opType.compare("explicit") == 0 ) {
      haveExplicitParabolicTerms = true;
      thermalDiffusionStatus.isExplicit = true;
    } else if(opType.compare("rkl") == 0 ) {
      haveRKLParabolicTerms = true;
      thermalDiffusionStatus.isRKL = true;
    } else {
      std::stringstream msg;
      msg  << "Unknown integration type for thermal diffusion: " << opType;
      IDEFIX_ERROR(msg);
    }
    this->thermalDiffusion.Init(input, grid, this);
  }

#if MHD == YES
  if(input.CheckEntry(std::string(Phys::prefix),"resistivity")>=0 ||
     input.CheckEntry(std::string(Phys::prefix),"ambipolar")>=0 ||
     input.CheckEntry(std::string(Phys::prefix),"hall")>=0 ) {
    //
    this->haveCurrent = true;

    if(input.CheckEntry(std::string(Phys::prefix),"resistivity")>=0) {
      std::string opType = input.Get<std::string>(std::string(Phys::prefix),"resistivity",0);
      if(opType.compare("explicit") == 0 ) {
        haveExplicitParabolicTerms = true;
        resistivityStatus.isExplicit = true;
        needExplicitCurrent = true;
      } else if(opType.compare("rkl") == 0 ) {
        haveRKLParabolicTerms = true;
        resistivityStatus.isRKL = true;
        needRKLCurrent = true;
      } else {
        std::stringstream msg;
        msg  << "Unknown integration type for resistivity: " << opType;
        IDEFIX_ERROR(msg);
      }
      if(input.Get<std::string>(std::string(Phys::prefix),"resistivity",1).compare("constant") == 0) {
        this->etaO = input.Get<real>(std::string(Phys::prefix),"resistivity",2);
        resistivityStatus.status = Constant;
      } else if(input.Get<std::string>(std::string(Phys::prefix),"resistivity",1).compare("userdef") == 0) {
        resistivityStatus.status = UserDefFunction;
      } else {
        IDEFIX_ERROR("Unknown resistivity definition in idefix.ini. "
                     "Can only be constant or userdef.");
      }
    }

    if(input.CheckEntry(std::string(Phys::prefix),"ambipolar")>=0) {
      std::string opType = input.Get<std::string>(std::string(Phys::prefix),"ambipolar",0);
      if(opType.compare("explicit") == 0 ) {
        haveExplicitParabolicTerms = true;
        ambipolarStatus.isExplicit = true;
        needExplicitCurrent = true;
      } else if(opType.compare("rkl") == 0 ) {
        haveRKLParabolicTerms = true;
        ambipolarStatus.isRKL = true;
        needRKLCurrent = true;
      } else {
        std::stringstream msg;
        msg  << "Unknown integration type for ambipolar: " << opType;
        IDEFIX_ERROR(msg);
      }
      if(input.Get<std::string>(std::string(Phys::prefix),"ambipolar",1).compare("constant") == 0) {
        this->xA = input.Get<real>(std::string(Phys::prefix),"ambipolar",2);
        ambipolarStatus.status = Constant;
      } else if(input.Get<std::string>(std::string(Phys::prefix),"ambipolar",1).compare("userdef") == 0) {
        ambipolarStatus.status = UserDefFunction;
      } else {
        IDEFIX_ERROR("Unknown ambipolar definition in idefix.ini. "
                     "Can only be constant or userdef.");
      }
    }

    if(input.CheckEntry(std::string(Phys::prefix),"hall")>=0) {
      std::string opType = input.Get<std::string>(std::string(Phys::prefix),"hall",0);
      if(opType.compare("explicit") == 0 ) {
        hallStatus.isExplicit = true;
        needExplicitCurrent = true;
      } else if(opType.compare("rkl") == 0 ) {
        IDEFIX_ERROR("RKL inegration is incompatible with Hall");
      } else {
        std::stringstream msg;
        msg  << "Unknown integration type for hall: " << opType;
        IDEFIX_ERROR(msg);
      }
      if(input.Get<std::string>(std::string(Phys::prefix),"hall",1).compare("constant") == 0) {
        this->xH = input.Get<real>(std::string(Phys::prefix),"hall",2);
        hallStatus.status = Constant;
      } else if(input.Get<std::string>(std::string(Phys::prefix),"hall",1).compare("userdef") == 0) {
        hallStatus.status = UserDefFunction;
      } else {
        IDEFIX_ERROR("Unknown Hall definition in idefix.ini. Can only be constant or userdef.");
      }
    }
  }
  #endif // MHD

  /////////////////////////////////////////
  //  ALLOCATION SECION ///////////////////
  /////////////////////////////////////////

  // We now allocate the fields required by the hydro solver
  Vc = IdefixArray4D<real>("FLUID_Vc", NVAR,
                           data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  Uc = IdefixArray4D<real>("FLUID_Uc", NVAR,
                           data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

  data->states["current"].PushArray(Uc, State::center, "FLUID_Uc");

  InvDt = IdefixArray3D<real>("FLUID_InvDt",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  cMax = IdefixArray3D<real>("FLUID_cMax",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  dMax = IdefixArray3D<real>("FLUID_dMax",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  FluxRiemann =  IdefixArray4D<real>("FLUID_FluxRiemann", NVAR,
                                     data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

    if constexpr(Phys::mhd) {
      Vs = IdefixArray4D<real>("FLUID_Vs", DIMENSIONS,
                  data->np_tot[KDIR]+KOFFSET, data->np_tot[JDIR]+JOFFSET, data->np_tot[IDIR]+IOFFSET);
      #ifdef EVOLVE_VECTOR_POTENTIAL
        #if DIMENSIONS == 1
          IDEFIX_ERROR("EVOLVE_VECTOR_POTENTIAL is not compatible with 1D MHD");
        #else
          Ve = IdefixArray4D<real>("FLUID_Ve", AX3e+1,
                data->np_tot[KDIR]+KOFFSET, data->np_tot[JDIR]+JOFFSET, data->np_tot[IDIR]+IOFFSET);

          data->states["current"].PushArray(Ve, State::center, "FLUID_Ve");
        #endif
      #else // EVOLVE_VECTOR_POTENTIAL
        data->states["current"].PushArray(Vs, State::center, "FLUID_Vs");
      #endif // EVOLVE_VECTOR_POTENTIAL
      
    }

  // Allocate sound speed array if needed
  if(this->haveIsoSoundSpeed == UserDefFunction) {
    this->isoSoundSpeedArray = IdefixArray3D<real>("FLUID_csIso",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  }

  if(this->haveCurrent) {
    // Allocate current (when hydro needs it)
    J = IdefixArray4D<real>("FLUID_J", 3,
                            data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  }

  // Allocate nonideal MHD effects array when a user-defined function is used
  if(this->resistivityStatus.status ==  UserDefFunction)
    etaOhmic = IdefixArray3D<real>("FLUID_etaOhmic",
                                    data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  if(this->ambipolarStatus.status == UserDefFunction)
    xAmbipolar = IdefixArray3D<real>("FLUID_xAmbipolar",
                                     data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  if(this->hallStatus.status == UserDefFunction)
    xHall = IdefixArray3D<real>("FLUID_xHall",
                                  data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

  // Fill the names of the fields
  for(int i = 0 ; i < NVAR ;  i++) {
    switch(i) {
      case RHO:
        VcName.push_back("RHO");
        data->vtk->RegisterVariable(Vc, "RHO", RHO);
        break;
      case VX1:
        VcName.push_back("VX1");
        data->vtk->RegisterVariable(Vc, "VX1", VX1);
        break;
      case VX2:
        VcName.push_back("VX2");
        data->vtk->RegisterVariable(Vc, "VX2", VX2);
        break;
      case VX3:
        VcName.push_back("VX3");
        data->vtk->RegisterVariable(Vc, "VX3", VX3);
        break;
      case BX1:
        VcName.push_back("BX1");
        data->vtk->RegisterVariable(Vc, "BX1", BX1);
        break;
      case BX2:
        VcName.push_back("BX2");
        data->vtk->RegisterVariable(Vc, "BX2", BX2);
        break;
      case BX3:
        VcName.push_back("BX3");
        data->vtk->RegisterVariable(Vc, "BX3", BX3);
        break;
#if HAVE_ENERGY
      case PRS:
        VcName.push_back("PRS");
        data->vtk->RegisterVariable(Vc, "PRS", PRS);
        break;
#endif
      default:
        VcName.push_back("Vc_"+std::to_string(i));
        data->vtk->RegisterVariable(Vc, "Vc_"+std::to_string(i), i);
    }
  }

  for(int i = 0 ; i < DIMENSIONS ; i++) {
    switch(i) {
      case 0:
        VsName.push_back("BX1s");
        break;
      case 1:
        VsName.push_back("BX2s");
        break;
      case 2:
        VsName.push_back("BX3s");
        break;
      default:
        VsName.push_back("Vs_"+std::to_string(i));
    }
  }

  #ifdef EVOLVE_VECTOR_POTENTIAL
    #if DIMENSIONS < 3
      VeName.push_back("AX3e");
    #else
      for(int i = 0 ; i < DIMENSIONS ; i++) {
      switch(i) {
        case 0:
          VeName.push_back("AX1e");
          break;
        case 1:
          VeName.push_back("AX2e");
          break;
        case 2:
          VeName.push_back("AX3e");
          break;
        default:
          VeName.push_back("Ve_"+std::to_string(i));
        }
      }
    #endif
  #endif


  //*******************************************
  //** Child object allocation section
  //*********************************************

    // Initialise Riemann Solver
  this->rSolver = std::unique_ptr<RiemannSolver<Phys>> (new RiemannSolver(input, this));

  if constexpr(Phys::mhd) {
    this->emf = std::unique_ptr<ConstrainedTransport<Phys>>(new ConstrainedTransport<Phys>(input, this));
  }


  // Do we have to take care of the axis?
  if(data->haveAxis) {
    this->myAxis = std::unique_ptr<Axis<Phys>>(new Axis<Phys>(grid, this));
    this->haveAxis = true;
  }

  // Initialise boundary conditions
  boundary = std::unique_ptr<Boundary<Phys>> (new Boundary<Phys>(input, grid, this));

  if(haveRKLParabolicTerms) {
    this->rkl = std::unique_ptr<RKLegendre<Phys>> (new RKLegendre<Phys>(input,this));
  }

  idfx::popRegion();
}


#include "addSourceTerms.hpp"
#include "calcRightHandSide.hpp"
#include "enroll.hpp"
#include "calcCurrent.hpp"
#include "coarsenFlow.hpp"
#include "convertConsToPrim.hpp"
#include "checkDivB.hpp"
#include "evolveStage.hpp"
#include "showConfig.hpp"
#endif // FLUID_FLUID_HPP_
