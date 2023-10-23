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
#include <memory>

#include "idefix.hpp"
#include "grid.hpp"
#include "fluid_defs.hpp"
#include "eos.hpp"
#include "thermalDiffusion.hpp"
#include "bragThermalDiffusion.hpp"
#include "selfGravity.hpp"
#include "vtk.hpp"
#include "dump.hpp"
#ifdef WITH_HDF5
#include "xdmf.hpp"
#endif

// forward class declaration
class DataBlock;
template<typename Phys>
class Boundary;

template<typename Phys>
class ConstrainedTransport;

template<typename Phys>
class RKLegendre;

template<typename Phys>
class RiemannSolver;

template<typename Phys>
class ShockFlattening;

class Viscosity;
class ThermalDiffusion;
class BragViscosity;
class BragThermalDiffusion;
class Drag;
class Tracer;


template<typename Phys>
class Fluid {
 public:
  Fluid( Grid &, Input&, DataBlock *, int n = 0);
  void ConvertConsToPrim();
  void ConvertPrimToCons();
  template <int> void CalcParabolicFlux(const real);
  template <int> void AddNonIdealMHDFlux(const real);
  template <int> void CalcRightHandSide(real, real );
  void CalcCurrent();
  void AddSourceTerms(real, real );
  void CoarsenFlow(IdefixArray4D<real>&);
  void CoarsenMagField(IdefixArray4D<real>&);
  real CheckDivB();
  void EvolveStage(const real, const real);
  void ResetStage();
  void ShowConfig();
  IdefixArray4D<real> GetFlux() {return this->FluxRiemann;}
  int CheckNan();

  // Our boundary conditions
  std::unique_ptr<Boundary<Phys>> boundary;

  // EOS
  std::unique_ptr<EquationOfState> eos;

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
  std::unique_ptr<ThermalDiffusion> thermalDiffusion;

  // Whether or not we have braginskii viscosity
  ParabolicModuleStatus bragViscosityStatus;

  // Whether or not we have braginskii thermal diffusion
  ParabolicModuleStatus bragThermalDiffusionStatus;

  // BragViscosity object
  std::unique_ptr<BragViscosity> bragViscosity;

  // Braginskii Thermal Diffusion object
  std::unique_ptr<BragThermalDiffusion> bragThermalDiffusion;

  // Drag object
  bool haveDrag{false};
  std::unique_ptr<Drag> drag;

  // Whether or not we have to treat the axis
  bool haveAxis{false};

  // Rotation vector
  bool haveRotation{false};
  real OmegaZ;

  bool haveShearingBox{false};
  // Shear rate for shearing box problems
  real sbS;
  // Box width for shearing box problems
  real sbLx;

  // Tracers treatment
  std::unique_ptr<Tracer> tracer;
  bool haveTracer{false};
  int nTracer{0};


  // Enroll user-defined boundary conditions (proxies for boundary class functions)
  template <typename T>
  void EnrollUserDefBoundary(T);
  template <typename T>
  void EnrollInternalBoundary(T);
  template <typename T>
  void EnrollFluxBoundary(T);

  void EnrollEmfBoundary(EmfBoundaryFunc);

  // Add some user source terms
  void EnrollUserSourceTerm(SrcTermFunc<Phys>);
  void EnrollUserSourceTerm(SrcTermFuncOld); // Deprecated

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

  DataBlock *data;

  // Data related to current instance of the Fluid object
  std::string prefix;
  int instanceNumber;

 private:
  friend class ConstrainedTransport<Phys>;
  friend class Fargo;
  friend class RKLegendre<Phys>;
  friend class Boundary<Phys>;
  friend class ShockFlattening<Phys>;
  friend class RiemannSolver<Phys>;
  friend class Viscosity;
  friend class ThermalDiffusion;
  friend class BragViscosity;
  friend class BragThermalDiffusion;
  friend class Drag;

  template <typename P>
  friend struct Fluid_AddSourceTermsFunctor;

  template <typename P, int dir>
  friend struct Fluid_CorrectFluxFunctor;

  template <typename P, int dir>
  friend struct Fluid_CalcRHSFunctor;

  template<typename P>
  friend struct ShockFlattening_FindShockFunctor;

  // Emf boundary conditions
  bool haveEmfBoundary{false};
  EmfBoundaryFunc emfBoundaryFunc{NULL};

  // User defined source term
  SrcTermFunc<Phys> userSourceTerm{NULL};
  SrcTermFuncOld    userSourceTermOld{NULL};
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

#include "physics.hpp"
#include "dataBlock.hpp"
#include "boundary.hpp"
#include "constrainedTransport.hpp"
#include "axis.hpp"
#include "rkl.hpp"
#include "riemannSolver.hpp"
#include "viscosity.hpp"
#include "bragViscosity.hpp"
#include "drag.hpp"
#include "checkNan.hpp"
#include "tracer.hpp"


template<typename Phys>
Fluid<Phys>::Fluid(Grid &grid, Input &input, DataBlock *datain, int n) {
  idfx::pushRegion("Fluid::Fluid");
  // Save the datablock to which we are attached from now on
  this->data = datain;

  // Create our own prefix
  prefix = std::string(Phys::prefix);

  // When dealing with dust, add the specie number
  if(Phys::prefix.compare("Dust") == 0) prefix += std::to_string(n);

  // Keep the instance # for later use
  instanceNumber = n;

  #if ORDER < 1 || ORDER > 4
     IDEFIX_ERROR("Reconstruction at chosen order is not implemented. Check your definitions file");
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

  // Passive tracers
  if(input.CheckEntry(std::string(Phys::prefix),"tracer")>=0) {
    this->haveTracer = true;
    // nTracer is initialised before instanciation of the Tracer object
    // Because we need to know the number of tracers to allocate Vc, Uc and Flux.
    this->nTracer = input.Get<int>(std::string(Phys::prefix),"tracer",0);
    if(this->nTracer < 1) {
      IDEFIX_ERROR("The number of passive tracers should be >= 1");
    }
  }

  // If we are not the primary hydro object, we copy the properties of the primary hydro object
  // so that we solve for consistant physics
  if(prefix.compare("Hydro") != 0) {
    this->haveSourceTerms = data->hydro->haveSourceTerms;
    this->haveRotation = data->hydro->haveRotation;
    this->OmegaZ = data->hydro->OmegaZ;
    this->haveShearingBox = data->hydro->haveShearingBox;
    this->sbS = data->hydro->sbS;
    this->sbLx = data->hydro->sbLx;
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
    if(input.Get<std::string>(std::string(Phys::prefix),"viscosity",1)
                              .compare("constant") == 0) {
      viscosityStatus.status = Constant;
    } else if(input.Get<std::string>(std::string(Phys::prefix),"viscosity",1)
                                      .compare("userdef") == 0) {
        viscosityStatus.status = UserDefFunction;
    } else {
        IDEFIX_ERROR("Unknown viscosity definition in idefix.ini. "
                     "Can only be constant or userdef.");
    }
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
    if(input.Get<std::string>(std::string(Phys::prefix),"TDiffusion",1).compare("constant") == 0) {
        thermalDiffusionStatus.status = Constant;
      } else if(input.Get<std::string>(std::string(Phys::prefix),
                                      "TDiffusion",1).compare("userdef") == 0) {
       thermalDiffusionStatus.status = UserDefFunction;

      } else {
        IDEFIX_ERROR("Unknown thermal diffusion definition in idefix.ini. "
                     "Can only be constant or userdef.");
      }
  }

  // Check whether braginskii viscosity is enabled, if so, init a braginskii viscosity object
  if(input.CheckEntry(std::string(Phys::prefix),"bragViscosity")>=0) {
    std::string opType = input.Get<std::string>(std::string(Phys::prefix),"bragViscosity",0);
    if(opType.compare("explicit") == 0 ) {
      haveExplicitParabolicTerms = true;
      bragViscosityStatus.isExplicit = true;
    } else if(opType.compare("rkl") == 0 ) {
      haveRKLParabolicTerms = true;
      bragViscosityStatus.isRKL = true;
    } else {
      std::stringstream msg;
      msg  << "Unknown integration type for braginskii viscosity: " << opType;
      IDEFIX_ERROR(msg);
    }
    this->bragViscosity = std::make_unique<BragViscosity>(input, grid, this);
  }

  // Check whether braginskii thermal diffusion is enabled,
  // if so, init a braginskii thermal diffusion object
  if(input.CheckEntry(std::string(Phys::prefix),"bragTDiffusion")>=0) {
    std::string opType = input.Get<std::string>(std::string(Phys::prefix),"bragTDiffusion",0);
    if(opType.compare("explicit") == 0 ) {
      haveExplicitParabolicTerms = true;
      bragThermalDiffusionStatus.isExplicit = true;
    } else if(opType.compare("rkl") == 0 ) {
      haveRKLParabolicTerms = true;
      bragThermalDiffusionStatus.isRKL = true;
    } else {
      std::stringstream msg;
      msg  << "Unknown integration type for braginskii thermal diffusion: " << opType;
      IDEFIX_ERROR(msg);
    }
    this->bragThermalDiffusion = std::make_unique<BragThermalDiffusion>(input, grid, this);
  }

  if(input.CheckEntry(std::string(Phys::prefix),"drag")>=0) {
    haveDrag = true;
  }

  if constexpr(Phys::mhd) {
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
        if(input.Get<std::string>(
            std::string(Phys::prefix),"resistivity",1).compare("constant") == 0) {
          this->etaO = input.Get<real>(std::string(Phys::prefix),"resistivity",2);
          resistivityStatus.status = Constant;
        } else if(input.Get<std::string>(
            std::string(Phys::prefix),"resistivity",1).compare("userdef") == 0) {
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
        if(input.Get<std::string>(
            std::string(Phys::prefix),"ambipolar",1).compare("constant") == 0) {
          this->xA = input.Get<real>(std::string(Phys::prefix),"ambipolar",2);
          ambipolarStatus.status = Constant;
        } else if(input.Get<std::string>(
                    std::string(Phys::prefix),"ambipolar",1).compare("userdef") == 0) {
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
        } else if(input.Get<std::string>(
                    std::string(Phys::prefix),"hall",1).compare("userdef") == 0) {
          hallStatus.status = UserDefFunction;
        } else {
          IDEFIX_ERROR("Unknown Hall definition in idefix.ini. Can only be constant or userdef.");
        }
      }
    }
  } // MHD

  /////////////////////////////////////////
  //  ALLOCATION SECION ///////////////////
  /////////////////////////////////////////

  // We now allocate the fields required by the hydro solver
  Vc = IdefixArray4D<real>(prefix+"_Vc", Phys::nvar+nTracer,
                           data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  Uc = IdefixArray4D<real>(prefix+"_Uc", Phys::nvar+nTracer,
                           data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

  data->states["current"].PushArray(Uc, State::center, prefix+"_Uc");

  InvDt = IdefixArray3D<real>(prefix+"_InvDt",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  cMax = IdefixArray3D<real>(prefix+"_cMax",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  dMax = IdefixArray3D<real>(prefix+"_dMax",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  FluxRiemann =  IdefixArray4D<real>(prefix+"_FluxRiemann", Phys::nvar+nTracer,
                                     data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

  if constexpr(Phys::mhd) {
    Vs = IdefixArray4D<real>(prefix+"_Vs", DIMENSIONS,
              data->np_tot[KDIR]+KOFFSET, data->np_tot[JDIR]+JOFFSET, data->np_tot[IDIR]+IOFFSET);
    #ifdef EVOLVE_VECTOR_POTENTIAL
      #if DIMENSIONS == 1
        IDEFIX_ERROR("EVOLVE_VECTOR_POTENTIAL is not compatible with 1D MHD");
      #else
        Ve = IdefixArray4D<real>(prefix+"_Ve", AX3e+1,
              data->np_tot[KDIR]+KOFFSET, data->np_tot[JDIR]+JOFFSET, data->np_tot[IDIR]+IOFFSET);

        data->states["current"].PushArray(Ve, State::center, prefix+"_Ve");
      #endif
    #else // EVOLVE_VECTOR_POTENTIAL
      data->states["current"].PushArray(Vs, State::center, prefix+"_Vs");
    #endif // EVOLVE_VECTOR_POTENTIAL
  }

  if(this->haveCurrent) {
    // Allocate current (when hydro needs it)
    J = IdefixArray4D<real>(prefix+"_J", 3,
                            data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  }

  // Allocate nonideal MHD effects array when a user-defined function is used
  if(this->resistivityStatus.status ==  UserDefFunction)
    etaOhmic = IdefixArray3D<real>(prefix+"_etaOhmic",
                                    data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  if(this->ambipolarStatus.status == UserDefFunction)
    xAmbipolar = IdefixArray3D<real>(prefix+"_xAmbipolar",
                                     data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  if(this->hallStatus.status == UserDefFunction)
    xHall = IdefixArray3D<real>(prefix+"_xHall",
                                  data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

  // Fill the names of the fields
  std::string outputPrefix("");
  // If we have hydro, the output prefix is "" for backward compatibility
  if(prefix.compare("Hydro") != 0) {
    outputPrefix = prefix;
    outputPrefix += "_";
  }

  for(int i = 0 ; i < Phys::nvar+nTracer ;  i++) {
    switch(i) {
      case RHO:
        VcName.push_back("RHO");
        data->dump->RegisterVariable(Vc, outputPrefix+"Vc-RHO", RHO);
        break;
      case VX1:
        VcName.push_back("VX1");
        data->dump->RegisterVariable(Vc, outputPrefix+"Vc-VX1", VX1, IDIR);
        break;
      case VX2:
        VcName.push_back("VX2");
        data->dump->RegisterVariable(Vc, outputPrefix+"Vc-VX2", VX2, JDIR);
        break;
      case VX3:
        VcName.push_back("VX3");
        data->dump->RegisterVariable(Vc, outputPrefix+"Vc-VX3", VX3, KDIR);
        break;
      case BX1:
        VcName.push_back("BX1");
        // never save cell-centered BX1 in dumps
        break;
      case BX2:
        VcName.push_back("BX2");
        #if DIMENSIONS < 2
          data->dump->RegisterVariable(Vc, outputPrefix+"Vc-BX2", BX2, JDIR);
        #endif
        break;
      case BX3:
        VcName.push_back("BX3");
        #if DIMENSIONS < 3
          data->dump->RegisterVariable(Vc, outputPrefix+"Vc-BX3", BX3, KDIR);
        #endif
        break;
      case PRS:
        VcName.push_back("PRS");
        data->dump->RegisterVariable(Vc, outputPrefix+"Vc-PRS", PRS);
        break;
      default:
        if(i>=Phys::nvar) {
          std::string tracerLabel = std::string("TR")+std::to_string(i-Phys::nvar); // ="TRn"
          VcName.push_back(tracerLabel);
          data->dump->RegisterVariable(Vc, outputPrefix+"Vc-"+tracerLabel, i);
        } else {
          VcName.push_back("Vc-"+std::to_string(i));
          data->vtk->RegisterVariable(Vc, outputPrefix+"Vc-"+std::to_string(i), i);
        }
    }
    data->vtk->RegisterVariable(Vc, outputPrefix+VcName[i], i);
    #ifdef WITH_HDF5
      data->xdmf->RegisterVariable(Vc, outputPrefix+VcName[i], i);
    #endif
  }

  if constexpr(Phys::mhd) {
      for(int i = 0 ; i < DIMENSIONS ; i++) {
        switch(i) {
          case 0:
            VsName.push_back("BX1s");
            data->dump->RegisterVariable(Vs, outputPrefix+"Vs-BX1s",
                                        BX1s, IDIR, DumpField::ArrayLocation::Face);
            break;
          case 1:
            VsName.push_back("BX2s");
            data->dump->RegisterVariable(Vs, outputPrefix+"Vs-BX2s",
                                        BX2s, JDIR, DumpField::ArrayLocation::Face);
            break;
          case 2:
            VsName.push_back("BX3s");
            data->dump->RegisterVariable(Vs, outputPrefix+"Vs-BX3s",
                                        BX3s, KDIR, DumpField::ArrayLocation::Face);
            break;
          default:
            VsName.push_back("Vs-"+std::to_string(i));
        }
      }

    #ifdef EVOLVE_VECTOR_POTENTIAL
      #if DIMENSIONS < 3
        VeName.push_back("AX3e");
        data->dump->RegisterVariable(Ve, outputPrefix+"Ve-AX3e",
                                      AX3e, KDIR, DumpField::ArrayLocation::Edge);
      #else
        for(int i = 0 ; i < DIMENSIONS ; i++) {
        switch(i) {
          case 0:
            VeName.push_back("AX1e");
            data->dump->RegisterVariable(Ve, outputPrefix+"Ve-AX1e",
                                        AX1e, IDIR, DumpField::ArrayLocation::Edge);
            break;
          case 1:
            VeName.push_back("AX2e");
            data->dump->RegisterVariable(Ve, outputPrefix+"Ve-AX2e",
                                        AX2e, JDIR, DumpField::ArrayLocation::Edge);
            break;
          case 2:
            VeName.push_back("AX3e");
            data->dump->RegisterVariable(Ve, outputPrefix+"Ve-AX3e",
                                        AX3e, KDIR, DumpField::ArrayLocation::Edge);
            break;
          default:
            VeName.push_back("Ve-"+std::to_string(i));
          }
        }
      #endif
    #endif
  }


  //*******************************************
  //** Child object allocation section
  //*********************************************

    // Initialise Riemann Solver
  this->rSolver = std::make_unique<RiemannSolver<Phys>>(input, this);

  if constexpr(Phys::mhd) {
    this->emf = std::make_unique<ConstrainedTransport<Phys>>(input, this);
  }

  // Initialise the EOS
  if constexpr(Phys::eos) {
    this->eos = std::make_unique<EquationOfState>(input, data, this->prefix);
  }

  // Initialise boundary conditions
  boundary = std::make_unique<Boundary<Phys>>(this);
  this->haveAxis = data->haveAxis;

  if(haveRKLParabolicTerms) {
    this->rkl = std::make_unique<RKLegendre<Phys>>(input,this);
  }

  // Thermal diffusion
  if(thermalDiffusionStatus.status != Disabled ) {
    this->thermalDiffusion = std::make_unique<ThermalDiffusion>(input, grid, this);
  }

  // Viscosity
  if(viscosityStatus.status != Disabled) {
    this->viscosity = std::make_unique<Viscosity>(input, grid, this);
  }

  // Braginskii Thermal diffusion
  if(bragThermalDiffusionStatus.status != Disabled ) {
    this->bragThermalDiffusion = std::make_unique<BragThermalDiffusion>(input, grid, this);
  }

  // Braginskii Viscosity
  if(bragViscosityStatus.status != Disabled) {
    this->bragViscosity = std::make_unique<BragViscosity>(input, grid, this);
  }


  // Drag force when needed
  if(haveDrag) {
    this->drag = std::make_unique<Drag>(input, this);
  }

  // Tracers when needed
  if(haveTracer) {
    this->tracer= std::make_unique<Tracer>(this, nTracer);
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
