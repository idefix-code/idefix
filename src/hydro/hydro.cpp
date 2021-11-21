// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
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

  // Initialise boundary conditions
  boundary.Init(input, grid, this);

  if(input.CheckEntry("Hydro","gamma")>0) {
    this->gamma = input.GetReal("Hydro","gamma",0);
    idfx::cout << "Hydro: adiabatic EOS with gamma=" << this->gamma <<std::endl;
  } else {
    this->gamma = 5.0/3.0;
#if HAVE_ENERGY
    idfx::cout << "Hydro: Warning! no gamma has been set in the input file, assuming gamma=5/3."
               << std::endl;
#endif
  }

  if(input.CheckEntry("Hydro","csiso")>0) {
    #if HAVE_ENERGY
      IDEFIX_ERROR("csiso is useless without the isothermal EOS enabled in definition.h");
    #else
      std::string isoString = input.GetString("Hydro","csiso",0);
      if(isoString.compare("constant") == 0) {
        this->haveIsoSoundSpeed = Constant;
        this->isoSoundSpeed = input.GetReal("Hydro","csiso",1);
      } else if(isoString.compare("userdef") == 0) {
        this->haveIsoSoundSpeed = UserDefFunction;
      } else {
        IDEFIX_ERROR("csiso admits only constant or userdef entries");
      }
    #endif
  } else {
#if HAVE_ENERGY
    // set the isothermal soundspeed, even though it will not be used
    this->isoSoundSpeed = -1.0;
#else
    IDEFIX_ERROR("You are using the ISOTHERMAL approximation "
                 "but have not set csiso in the ini file.");
#endif
  }

  // read Solver from input file
  std::string solverString = input.GetString("Hydro","solver",0);

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
    this->OmegaZ = input.GetReal("Hydro","rotation",0);

    idfx::cout << "Hydro: Rotation enabled with Omega=" << this->OmegaZ << std::endl;
  }

  // Check whether we have shearing box
  int shearingbox = input.CheckEntry("Hydro","shearingBox");

  if(shearingbox>=0 ) {
    this->haveShearingBox = true;
    this->haveSourceTerms = true;
    if(shearingbox != 1) {
      IDEFIX_ERROR("Shearing box needs a scalar value for the shear rate in idefix.ini");
    }

    this->sbS = input.GetReal("Hydro","shearingBox",0);
    // Get box size
    this->sbLx = grid.xend[IDIR] - grid.xbeg[IDIR];

    idfx::cout << "Hydro: ShearingBox enabled with Shear rate= " << this->sbS
               << " and Lx= " << sbLx << std::endl;
  }

  // Gravitational potential (deprecated, use [Gravity] block in your input file)
  if(input.CheckEntry("Hydro","gravPotential")>=0) {
    IDEFIX_DEPRECATED("gravPotential in [Hydro] block is deprecated in the input file. "
                      "Use a [Gravity] block instead.");
    std::string potentialString = input.GetString("Hydro","gravPotential",0);
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
    std::string potentialString = input.GetString("Hydro","bodyForce",0);
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


  // Check whether viscosity is enabled, if so, construct a viscosity object
  if(input.CheckEntry("Hydro","viscosity")>=0) {
    std::string opType = input.GetString("Hydro","viscosity",0);
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

#if MHD == YES
  if(input.CheckEntry("Hydro","resistivity")>=0 ||
     input.CheckEntry("Hydro","ambipolar")>=0 ||
     input.CheckEntry("Hydro","hall")>=0 ) {
    //
    this->haveCurrent = true;

    if(input.CheckEntry("Hydro","resistivity")>=0) {
      std::string opType = input.GetString("Hydro","resistivity",0);
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
      if(input.GetString("Hydro","resistivity",1).compare("constant") == 0) {
        idfx::cout << "Hydro: Enabling Ohmic resistivity with constant diffusivity." << std::endl;
        this->etaO = input.GetReal("Hydro","resistivity",2);
        resistivityStatus.status = Constant;
      } else if(input.GetString("Hydro","resistivity",1).compare("userdef") == 0) {
        idfx::cout << "Hydro: Enabling Ohmic resistivity with user-defined diffusivity function."
                   << std::endl;
        resistivityStatus.status = UserDefFunction;
      } else {
        IDEFIX_ERROR("Unknown resistivity definition in idefix.ini. "
                     "Can only be constant or userdef.");
      }
    }

    if(input.CheckEntry("Hydro","ambipolar")>=0) {
      std::string opType = input.GetString("Hydro","ambipolar",0);
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
      if(input.GetString("Hydro","ambipolar",1).compare("constant") == 0) {
        idfx::cout << "Hydro: Enabling ambipolar diffusion with constant diffusivity."
                   << std::endl;
        this->xA = input.GetReal("Hydro","ambipolar",2);
        ambipolarStatus.status = Constant;
      } else if(input.GetString("Hydro","ambipolar",1).compare("userdef") == 0) {
        idfx::cout << "Hydro: Enabling ambipolar diffusion with user-defined diffusivity function."
                   << std::endl;
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
      #if EMF_AVERAGE != ARITHMETIC
        IDEFIX_ERROR("the Hall effect module is demonstrated stable only when using "
                    "EMF_AVERAGE=ARITHMETIC");
      #endif
      std::string opType = input.GetString("Hydro","hall",0);
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
      if(input.GetString("Hydro","hall",1).compare("constant") == 0) {
        idfx::cout << "Hydro: Enabling Hall effect with constant diffusivity." << std::endl;
        this->xH = input.GetReal("Hydro","hall",2);
        hallStatus.status = Constant;
      } else if(input.GetString("Hydro","hall",1).compare("userdef") == 0) {
        idfx::cout << "Hydro: Enabling Hall effect with user-defined diffusivity function."
                   << std::endl;
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
  Uc0 = IdefixArray4D<real>("Hydro_Uc0", NVAR,
                           data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

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

    Vs0 = IdefixArray4D<real>("Hydro_Vs0", DIMENSIONS,
                data->np_tot[KDIR]+KOFFSET, data->np_tot[JDIR]+JOFFSET, data->np_tot[IDIR]+IOFFSET);
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


  idfx::popRegion();
}

void Hydro::EnrollIsoSoundSpeed(IsoSoundSpeedFunc myFunc) {
  if(this->haveIsoSoundSpeed != UserDefFunction) {
    IDEFIX_ERROR("Isothermal sound speed enrollment requires Hydro/csiso "
                 " to be set to userdef in .ini file");
  }
  #if HAVE_ENERGY
    IDEFIX_ERROR("Isothermal sound speed enrollment requires ISOTHERMAL to be defined in"
                 "definitions.hpp");
  #endif
  this->isoSoundSpeedFunc = myFunc;
  idfx::cout << "Hydro: User-defined isothermal sound speed has been enrolled" << std::endl;
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
  this->emfBoundaryFunc = myFunc;
  this->haveEmfBoundary = true;
  idfx::cout << "Hydro: User-defined EMF boundary condition has been enrolled" << std::endl;
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
  idfx::cout << "Hydro: User-defined source term has been enrolled" << std::endl;
}

void Hydro::EnrollOhmicDiffusivity(DiffusivityFunc myFunc) {
  if(this->resistivityStatus.status < UserDefFunction) {
    IDEFIX_ERROR("Ohmic diffusivity enrollment requires Hydro/Resistivity "
                 "to be set to userdef in .ini file");
  }
  this->ohmicDiffusivityFunc = myFunc;
  idfx::cout << "Hydro: User-defined ohmic diffusivity has been enrolled" << std::endl;
}

void Hydro::EnrollAmbipolarDiffusivity(DiffusivityFunc myFunc) {
  if(this->ambipolarStatus.status < UserDefFunction) {
    IDEFIX_ERROR("Ambipolar diffusivity enrollment requires Hydro/Ambipolar "
                 "to be set to userdef in .ini file");
  }
  this->ambipolarDiffusivityFunc = myFunc;
  idfx::cout << "Hydro: User-defined ambipolar diffusivity has been enrolled" << std::endl;
}

void Hydro::EnrollHallDiffusivity(DiffusivityFunc myFunc) {
  if(this->hallStatus.status < UserDefFunction) {
    IDEFIX_ERROR("Hall diffusivity enrollment requires Hydro/Hall "
                 "to be set to userdef in .ini file");
  }
  this->hallDiffusivityFunc = myFunc;
  idfx::cout << "Hydro: User-defined Hall diffusivity has been enrolled" << std::endl;
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
