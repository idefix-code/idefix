// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>

#include "idefix.hpp"
#include "hydro.hpp"
#include "dataBlock.hpp"


void Hydro::Init(Input &input, Grid &grid, DataBlock *datain) {
  idfx::pushRegion("Hydro::Init");
  // Save the datablock to which we are attached from now on
  this->data = datain;

  #if ORDER < 1 || ORDER > 4
     IDEFIX_ERROR("Reconstruction at chosen order is not implemented. Check your definitions file");
  #endif


  #if HAVE_ENERGY
    this->gamma = input.GetOrSet<real>("Hydro","gamma",0, 5.0/3.0);
  #endif

  #ifdef ISOTHERMAL
    std::string isoString = input.Get<std::string>("Hydro","csiso",0);
    if(isoString.compare("constant") == 0) {
      this->haveIsoSoundSpeed = Constant;
      this->isoSoundSpeed = input.Get<real>("Hydro","csiso",1);
    } else if(isoString.compare("userdef") == 0) {
      this->haveIsoSoundSpeed = UserDefFunction;
    } else {
      IDEFIX_ERROR("csiso admits only constant or userdef entries");
    }
  #else
    // set the isothermal soundspeed, even though it will not be used
    this->isoSoundSpeed = -1.0;
  #endif

  // read Solver from input file
  std::string solverString = input.Get<std::string>("Hydro","solver",0);

  if (solverString.compare("tvdlf") == 0) {
    mySolver = TVDLF;
  } else if (solverString.compare("hll") == 0) {
    mySolver = HLL;
#if MHD == YES
  } else if (solverString.compare("hlld") == 0) {
    mySolver = HLLD;
#else
  } else if (solverString.compare("hllc") == 0) {
    mySolver = HLLC;
#endif
  } else if (solverString.compare("roe") == 0) {
    mySolver = ROE;
  } else {
    std::stringstream msg;
#if MHD == YES
    msg << "Unknown MHD solver type " << solverString;
#else
    msg << "Unknown HD solver type " << solverString;
#endif
    IDEFIX_ERROR(msg);
  }

  // Shock flattening
  this->haveShockFlattening = input.CheckEntry("Hydro","shockFlattening")>=0;

  // Source terms (always activated when non-cartesian geometry because of curvature source terms)
#if GEOMETRY == CARTESIAN
  this->haveSourceTerms = false;
#else
  this->haveSourceTerms = true;
#endif

  // Check whether we have rotation
  int rotation = input.CheckEntry("Hydro","rotation");

  if(rotation>=0 ) {
    this->haveSourceTerms = true;
    this->haveRotation = true;
    if(rotation >1 ) IDEFIX_ERROR("Rotation needs a a single rotation velocity (Omega_Z) \
                                   in idefix.ini");
    this->OmegaZ = input.Get<real>("Hydro","rotation",0);
  }

  // Check whether we have shearing box
  int shearingbox = input.CheckEntry("Hydro","shearingBox");

  if(shearingbox>=0 ) {
    this->haveShearingBox = true;
    this->haveSourceTerms = true;
    if(shearingbox != 1) {
      IDEFIX_ERROR("Shearing box needs a scalar value for the shear rate in idefix.ini");
    }

    this->sbS = input.Get<real>("Hydro","shearingBox",0);
    // Get box size
    this->sbLx = grid.xend[IDIR] - grid.xbeg[IDIR];
  }

  // Gravitational potential (deprecated, use [Gravity] block in your input file)
  if(input.CheckEntry("Hydro","gravPotential")>=0) {
    IDEFIX_DEPRECATED("gravPotential in [Hydro] block is deprecated in the input file. "
                      "Use a [Gravity] block instead.");
    std::string potentialString = input.Get<std::string>("Hydro","gravPotential",0);
    if(potentialString.compare("userdef") == 0) {
      data->gravity.haveUserDefPotential = true;
      data->gravity.havePotential = true;
      data->gravity.Init(input,data);
    } else {
      IDEFIX_ERROR("Unknown type of gravitational potential in idefix.ini. "
                   "Only userdef is implemented");
    }
  }

  // Body Force (deprecated, use [Gravity] block in your input file)
  if(input.CheckEntry("Hydro","bodyForce")>=0) {
    IDEFIX_DEPRECATED("bodyForce in [Hydro] block is deprecated in the input file. "
                      "Use a [Gravity] block instead.");
    std::string potentialString = input.Get<std::string>("Hydro","bodyForce",0);
    if(potentialString.compare("userdef") == 0) {
      data->gravity.haveBodyForce = true;
      data->gravity.Init(input,data);
    } else {
      IDEFIX_ERROR("Unknown type of body force in idefix.ini. "
                   "Only userdef is implemented");
    }
  }

  // Do we use fargo? (deprecated, use a full fargo block in idefix.ini,
  // instead of including it in hydo)
  if(input.CheckEntry("Hydro","fargo")>=0) {
    IDEFIX_DEPRECATED("Fargo in [Hydro] block is deprecated in the input file. "
                      "Use a [Fargo] block instead.");
    data->haveFargo = true;
    data->fargo.Init(input, this->data);
  }

  ///////////////////////
  // Parabolic terms
  ///////////////////////


  // Check whether viscosity is enabled, if so, init a viscosity object
  if(input.CheckEntry("Hydro","viscosity")>=0) {
    std::string opType = input.Get<std::string>("Hydro","viscosity",0);
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
    this->viscosity.Init(input, grid, this);
  }

  // Check whether thermal diffusion is enabled, if so, init a thermal diffusion object
  if(input.CheckEntry("Hydro","TDiffusion")>=0) {
    std::string opType = input.Get<std::string>("Hydro","TDiffusion",0);
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
  if(input.CheckEntry("Hydro","resistivity")>=0 ||
     input.CheckEntry("Hydro","ambipolar")>=0 ||
     input.CheckEntry("Hydro","hall")>=0 ) {
    //
    this->haveCurrent = true;

    if(input.CheckEntry("Hydro","resistivity")>=0) {
      std::string opType = input.Get<std::string>("Hydro","resistivity",0);
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
      if(input.Get<std::string>("Hydro","resistivity",1).compare("constant") == 0) {
        this->etaO = input.Get<real>("Hydro","resistivity",2);
        resistivityStatus.status = Constant;
      } else if(input.Get<std::string>("Hydro","resistivity",1).compare("userdef") == 0) {
        resistivityStatus.status = UserDefFunction;
      } else {
        IDEFIX_ERROR("Unknown resistivity definition in idefix.ini. "
                     "Can only be constant or userdef.");
      }
    }

    if(input.CheckEntry("Hydro","ambipolar")>=0) {
      std::string opType = input.Get<std::string>("Hydro","ambipolar",0);
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
      if(input.Get<std::string>("Hydro","ambipolar",1).compare("constant") == 0) {
        this->xA = input.Get<real>("Hydro","ambipolar",2);
        ambipolarStatus.status = Constant;
      } else if(input.Get<std::string>("Hydro","ambipolar",1).compare("userdef") == 0) {
        ambipolarStatus.status = UserDefFunction;
      } else {
        IDEFIX_ERROR("Unknown ambipolar definition in idefix.ini. "
                     "Can only be constant or userdef.");
      }
    }

    if(input.CheckEntry("Hydro","hall")>=0) {
      // Check consistency
      if(mySolver != HLL )
        IDEFIX_ERROR("Hall effect is only compatible with HLL Riemann solver.");
      std::string opType = input.Get<std::string>("Hydro","hall",0);
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
      if(input.Get<std::string>("Hydro","hall",1).compare("constant") == 0) {
        this->xH = input.Get<real>("Hydro","hall",2);
        hallStatus.status = Constant;
      } else if(input.Get<std::string>("Hydro","hall",1).compare("userdef") == 0) {
        hallStatus.status = UserDefFunction;
      } else {
        IDEFIX_ERROR("Unknown Hall definition in idefix.ini. Can only be constant or userdef.");
      }
    }
  }
  #endif // MHD

  // Do we have to take care of the axis?
  if(data->haveAxis) {
    this->myAxis.Init(grid, this);
    this->haveAxis = true;
  }
  /////////////////////////////////////////
  //  ALLOCATION SECION ///////////////////
  /////////////////////////////////////////

  // We now allocate the fields required by the hydro solver
  Vc = IdefixArray4D<real>("Hydro_Vc", NVAR,
                           data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  Uc = IdefixArray4D<real>("Hydro_Uc", NVAR,
                           data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

  data->states["current"].PushArray(Uc, State::center, "Hydro_Uc");

  InvDt = IdefixArray3D<real>("Hydro_InvDt",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  cMax = IdefixArray3D<real>("Hydro_cMax",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  dMax = IdefixArray3D<real>("Hydro_dMax",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  FluxRiemann =  IdefixArray4D<real>("Hydro_FluxRiemann", NVAR,
                                     data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

  #if MHD == YES
    Vs = IdefixArray4D<real>("Hydro_Vs", DIMENSIONS,
                data->np_tot[KDIR]+KOFFSET, data->np_tot[JDIR]+JOFFSET, data->np_tot[IDIR]+IOFFSET);
    #ifdef EVOLVE_VECTOR_POTENTIAL
      #if DIMENSIONS == 1
        IDEFIX_ERROR("EVOLVE_VECTOR_POTENTIAL is not compatible with 1D MHD");
      #else
        Ve = IdefixArray4D<real>("Hydro_Ve", AX3e+1,
              data->np_tot[KDIR]+KOFFSET, data->np_tot[JDIR]+JOFFSET, data->np_tot[IDIR]+IOFFSET);

        data->states["current"].PushArray(Ve, State::center, "Hydro_Ve");
      #endif
    #else // EVOLVE_VECTOR_POTENTIAL
      data->states["current"].PushArray(Vs, State::center, "Hydro_Vs");
    #endif // EVOLVE_VECTOR_POTENTIAL
    this->emf.Init(input, this);
  #endif

  // Allocate sound speed array if needed
  if(this->haveIsoSoundSpeed == UserDefFunction) {
    this->isoSoundSpeedArray = IdefixArray3D<real>("Hydro_csIso",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  }

  if(this->haveCurrent) {
    // Allocate current (when hydro needs it)
    J = IdefixArray4D<real>("Hydro_J", 3,
                            data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  }

  // Allocate nonideal MHD effects array when a user-defined function is used
  if(this->resistivityStatus.status ==  UserDefFunction)
    etaOhmic = IdefixArray3D<real>("Hydro_etaOhmic",
                                    data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  if(this->ambipolarStatus.status == UserDefFunction)
    xAmbipolar = IdefixArray3D<real>("Hydro_xAmbipolar",
                                     data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  if(this->hallStatus.status == UserDefFunction)
    xHall = IdefixArray3D<real>("Hydro_xHall",
                                  data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

  // Fill the names of the fields
  for(int i = 0 ; i < NVAR ;  i++) {
    switch(i) {
      case RHO:
        VcName.push_back("RHO");
        break;
      case VX1:
        VcName.push_back("VX1");
        break;
      case VX2:
        VcName.push_back("VX2");
        break;
      case VX3:
        VcName.push_back("VX3");
        break;
      case BX1:
        VcName.push_back("BX1");
        break;
      case BX2:
        VcName.push_back("BX2");
        break;
      case BX3:
        VcName.push_back("BX3");
        break;
#if HAVE_ENERGY
      case PRS:
        VcName.push_back("PRS");
        break;
#endif
      default:
        VcName.push_back("Vc_"+std::to_string(i));
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

  // Initialise boundary conditions
  boundary.Init(input, grid, this);

  // Init shock flattening
  if(haveShockFlattening) {
    this->shockFlattening = ShockFlattening(this,input.Get<real>("Hydro","shockFlattening",0));
  }

  idfx::popRegion();
}

void Hydro::EnrollIsoSoundSpeed(IsoSoundSpeedFunc myFunc) {
  if(this->haveIsoSoundSpeed != UserDefFunction) {
    IDEFIX_WARNING("Isothermal sound speed enrollment requires Hydro/csiso "
                 " to be set to userdef in .ini file");
  }
  #if HAVE_ENERGY
    IDEFIX_ERROR("Isothermal sound speed enrollment requires ISOTHERMAL to be defined in"
                 "definitions.hpp");
  #endif
  this->isoSoundSpeedFunc = myFunc;
}

void Hydro::EnrollUserDefBoundary(UserDefBoundaryFunc myFunc) {
  // This is a proxy for userdef enrollment
  boundary.EnrollUserDefBoundary(myFunc);
}

void Hydro::EnrollFluxBoundary(UserDefBoundaryFunc myFunc) {
  // This is a proxy for userdef enrollment
  boundary.EnrollFluxBoundary(myFunc);
}


void Hydro::EnrollInternalBoundary(InternalBoundaryFunc myFunc) {
  // This is a proxy for userdef enrollment
  boundary.EnrollInternalBoundary(myFunc);
}

void Hydro::EnrollEmfBoundary(EmfBoundaryFunc myFunc) {
  #if MHD == NO
    IDEFIX_ERROR("This function can only be used with the MHD solver.");
  #endif
  this->emfBoundaryFunc = myFunc;
  this->haveEmfBoundary = true;
}

// Deprecated function
void Hydro::EnrollGravPotential(GravPotentialFunc myFunc) {
  IDEFIX_DEPRECATED("Calling Hydro::EnrollGravPotential is deprecated."
                    "Use Gravity::EnrollPotential");
  data->gravity.EnrollPotential(myFunc);
}

// Deprecated function
void Hydro::EnrollBodyForce(BodyForceFunc myFunc) {
  IDEFIX_DEPRECATED("Calling Hydro::EnrollBodyForce is deprecated."
                    "Use Gravity::EnrollBodyForce");
  data->gravity.EnrollBodyForce(myFunc);
}

void Hydro::EnrollUserSourceTerm(SrcTermFunc myFunc) {
  this->userSourceTerm = myFunc;
  this->haveUserSourceTerm = true;
  this->haveSourceTerms = true;
}

void Hydro::EnrollOhmicDiffusivity(DiffusivityFunc myFunc) {
  #if MHD == NO
    IDEFIX_ERROR("This function can only be used with the MHD solver.");
  #endif
  if(this->resistivityStatus.status < UserDefFunction) {
    IDEFIX_WARNING("Ohmic diffusivity enrollment requires Hydro/Resistivity "
                 "to be set to userdef in .ini file");
  }
  this->ohmicDiffusivityFunc = myFunc;
}

void Hydro::EnrollAmbipolarDiffusivity(DiffusivityFunc myFunc) {
  #if MHD == NO
    IDEFIX_ERROR("This function can only be used with the MHD solver.");
  #endif
  if(this->ambipolarStatus.status < UserDefFunction) {
    IDEFIX_WARNING("Ambipolar diffusivity enrollment requires Hydro/Ambipolar "
                 "to be set to userdef in .ini file");
  }
  this->ambipolarDiffusivityFunc = myFunc;
}

void Hydro::EnrollHallDiffusivity(DiffusivityFunc myFunc) {
  #if MHD == NO
    IDEFIX_ERROR("This function can only be used with the MHD solver.");
  #endif
  if(this->hallStatus.status < UserDefFunction) {
    IDEFIX_WARNING("Hall diffusivity enrollment requires Hydro/Hall "
                 "to be set to userdef in .ini file");
  }
  this->hallDiffusivityFunc = myFunc;
}

Hydro::Hydro() {
  // Do nothing !
}

real Hydro::GetGamma() {
  return(this->gamma);
}


void Hydro::ResetStage() {
  // Reset variables required at the beginning of each stage
  // (essentially linked to timestep evaluation)
  idfx::pushRegion("Hydro::ResetStage");

  IdefixArray3D<real> InvDt=this->InvDt;

  idefix_for("HydroResetStage",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      InvDt(k,j,i) = ZERO_F;
  });

  idfx::popRegion();
}

void Hydro::ShowConfig() {
  idfx::cout << "Hydro: ";
  #if MHD == YES
    idfx::cout << "solving MHD equations." << std::endl;
    #ifdef EVOLVE_VECTOR_POTENTIAL
      idfx::cout << "Hydro: Using EXPERIMENTAL vector potential formulation for MHD." << std::endl;
    #endif
  #else
    idfx::cout << "solving HD equations." << std::endl;
  #endif
  idfx::cout << "Hydro: Reconstruction: ";
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
    idfx::cout << "Hydro: EOS: ideal with gamma=" << this->gamma << std::endl;
  #endif
  #ifdef ISOTHERMAL
    if(haveIsoSoundSpeed == Constant) {
      idfx::cout << "Hydro: EOS: isothermal with cs=" << isoSoundSpeed << "." << std::endl;
    } else if(haveIsoSoundSpeed == UserDefFunction) {
      idfx::cout << "Hydro: EOS: isothermal with user-defined cs function." << std::endl;
      if(!isoSoundSpeedFunc) {
        IDEFIX_ERROR("No user-defined isothermal sound speed function has been enrolled.");
      }
    }
  #endif// ISOTHERMAL
  idfx::cout << "Hydro: Riemann solver: ";
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
    idfx::cout << "Hydro: Rotation ENABLED with Omega=" << this->OmegaZ << std::endl;
  }
  if(haveShearingBox) {
    idfx::cout << "Hydro: ShearingBox ENABLED with shear rate= " << this->sbS
               << " and Lx= " << sbLx << std::endl;
  }
  if(haveShockFlattening) {
    idfx::cout << "Hydro: Shock Flattening ENABLED." << std::endl;
  }


  if(resistivityStatus.status != Disabled) {
    if(resistivityStatus.status == Constant) {
      idfx::cout << "Hydro: Ohmic resistivity ENABLED with constant resistivity eta="
                 << etaO << std::endl;
    } else if(resistivityStatus.status == UserDefFunction) {
      idfx::cout << "Hydro: Ohmic resistivity ENABLED with user-defined resistivity function."
                 << std::endl;
      if(!ohmicDiffusivityFunc) {
        IDEFIX_ERROR("No user-defined Ihmic resistivity function has been enrolled.");
      }
    } else {
      IDEFIX_ERROR("Unknown Ohmic resistivity mode");
    }
    if(resistivityStatus.isExplicit) {
      idfx::cout << "Hydro: Ohmic resistivity uses an explicit time integration." << std::endl;
    } else if(resistivityStatus.isRKL) {
      idfx::cout << "Hydro: Ohmic resistivity uses a Runge-Kutta-Legendre time integration."
                 << std::endl;
    } else {
      IDEFIX_ERROR("Unknown time integrator for Ohmic resistivity");
    }
  }

  if(ambipolarStatus.status != Disabled) {
    if(ambipolarStatus.status == Constant) {
      idfx::cout << "Hydro: Ambipolar diffusion ENABLED with constant diffusivity xA="
                 << xA << std::endl;
    } else if(ambipolarStatus.status == UserDefFunction) {
      idfx::cout << "Hydro: Ambipolar diffusion ENABLED with user-defined diffusivity function."
                 << std::endl;
      if(!ambipolarDiffusivityFunc) {
        IDEFIX_ERROR("No user-defined ambipolar diffusion function has been enrolled.");
      }
    } else {
      IDEFIX_ERROR("Unknown Ambipolar diffusion mode");
    }
    if(ambipolarStatus.isExplicit) {
      idfx::cout << "Hydro: Ambipolar diffusion uses an explicit time integration." << std::endl;
    } else if(ambipolarStatus.isRKL) {
      idfx::cout << "Hydro: Ambipolar diffusion uses a Runge-Kutta-Legendre time integration."
                 << std::endl;
    } else {
      IDEFIX_ERROR("Unknown time integrator for ambipolar diffusion");
    }
  }

  if(hallStatus.status != Disabled) {
    if(hallStatus.status == Constant) {
      idfx::cout << "Hydro: Hall effect ENABLED with constant diffusivity xH="
                 << xH << std::endl;
    } else if(hallStatus.status == UserDefFunction) {
      idfx::cout << "Hydro: Hall effect ENABLED with user-defined diffusivity function."
                 << std::endl;
      if(!hallDiffusivityFunc) {
        IDEFIX_ERROR("No user-defined Hall diffusivity function has been enrolled.");
      }
    } else {
      IDEFIX_ERROR("Unknown Hall effect mode");
    }
    if(hallStatus.isExplicit) {
      idfx::cout << "Hydro: Hall effect uses an explicit time integration." << std::endl;
    }  else {
      IDEFIX_ERROR("Unknown time integrator for Hall effect");
    }
  }

  if(emfBoundaryFunc) {
    idfx::cout << "Hydro: user-defined EMF boundaries ENABLED." << std::endl;
  }
  if(userSourceTerm) {
    idfx::cout << "Hydro: user-defined source terms ENABLED." << std::endl;
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
