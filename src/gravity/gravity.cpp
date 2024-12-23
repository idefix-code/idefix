// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>

#include "gravity.hpp"
#include "planetarySystem.hpp"
#include "dataBlock.hpp"
#include "output.hpp"
#include "input.hpp"

Gravity::Gravity(Input &input, DataBlock *datain) {
  idfx::pushRegion("Gravity::Gravity");
  this->data = datain;

  // Gravitational constant G
  this->gravCst = input.GetOrSet<real>("Gravity","gravCst",0, 1.0);
  // Gravitational potential
  int nPotential = input.CheckEntry("Gravity","potential");
  if(nPotential >=0) {
    this->havePotential = true;
    for(int i = 0 ; i < nPotential ; i++) {
      std::string potentialString = input.Get<std::string>("Gravity","potential",i);
      if(potentialString.compare("userdef") == 0) {
        this->haveUserDefPotential = true;
      } else if (potentialString.compare("central") == 0) {
        this->haveCentralMassPotential = true;
        this->centralMass = input.GetOrSet<real>("Gravity","Mcentral",0, 1.0);
        data->dump->RegisterVariable(&this->centralMass, "centralMass",1);
      } else if (potentialString.compare("selfgravity") == 0) {
        this->haveSelfGravityPotential = true;
      } else if (potentialString.compare("planet") == 0) {
        this->havePlanetsPotential = true;
        if(!datain->haveplanetarySystem) {
          IDEFIX_ERROR("You need to define a [Planet] block if"
                        "you want to work with embedded planets");
        }
      } else {
        IDEFIX_ERROR("Unknown type of gravitational potential in idefix.ini. ");
      }
    }
  }

  // Automatically enables gravity if a planetary system was initialised.
  if(datain->haveplanetarySystem) {
    this->havePlanetsPotential = true;
    this->havePotential = true;
  }

  // Body Force
  if(input.CheckEntry("Gravity","bodyForce")>=0) {
    std::string potentialString = input.Get<std::string>("Gravity","bodyForce",0);
    if(potentialString.compare("userdef") == 0) {
      this->haveBodyForce = true;
    } else {
      IDEFIX_ERROR("Unknown type of body force in idefix.ini. "
                   "Only userdef is implemented");
    }
  }



  // Allocate required arrays
  if(havePotential && (!haveInitialisedPotential)) {
    phiP = IdefixArray3D<real>("Gravity_PhiP",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
    haveInitialisedPotential = true;
  }
  if(haveBodyForce && !haveInitialisedBodyForce) {
    bodyForceVector = IdefixArray4D<real>("Gravity_bodyForce", COMPONENTS,
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
    haveInitialisedBodyForce = true;
  }

  // TODO(mauxionj): the following should consider the possibility to enable SelfGravity
  // even tought it is not explicitly asked by the user, to handle dependencies of other modules
  // (e.g. planet module ?). Those dependencies are to be handled in the previous conditional block
  // such that from this point, if SelfGravity is needed, haveSelfGravityPotential is true and
  // we build the corresponding object if not done yet.

  // Check SelfGravity object
  if(haveSelfGravityPotential) {
    selfGravity.Init(input, this->data);
    haveInitialisedSelfGravity = true;
  }

  this->skipGravity = input.GetOrSet<int>("Gravity","skip",0,1);
  if(skipGravity<1) {
    IDEFIX_ERROR("[Gravity]:skip should be a strictly positive integer");
  }
  idfx::popRegion();
}

void Gravity::ShowConfig() {
  if(data->haveGravity) {
    idfx::cout << "Gravity: ENABLED." << std::endl;
    idfx::cout << "Gravity: G=" << gravCst << "." << std::endl;
    if(haveUserDefPotential) {
      idfx::cout << "Gravity: User-defined gravitational potential ENABLED." << std::endl;
      if(!gravPotentialFunc) {
        IDEFIX_ERROR("No user-defined gravitational potential has been enrolled.");
      }
    }
    if(haveCentralMassPotential) {
      idfx::cout << "Gravity: central mass gravitational potential ENABLED with M="
                  << this->centralMass << std::endl;
    }
    if(haveSelfGravityPotential) {
      idfx::cout << "Gravity: self-gravity ENABLED." << std::endl;
      selfGravity.ShowConfig();
    }
    if(havePlanetsPotential) {
      idfx::cout << "Gravity: planet(s) potential ENABLED." << std::endl;
    }
    if(haveBodyForce) {
      idfx::cout << "Gravity: user-defined body force ENABLED." << std::endl;
      if(!bodyForceFunc) {
        IDEFIX_ERROR("No user-defined body force has been enrolled.");
      }
    }
    if(skipGravity>1) {
      idfx::cout << "Gravity: gravity field will be updated every " << skipGravity
                 << " cycles." << std::endl;
    }
  }
}
// This function compute the gravitational field, using both body force and potential
void Gravity::ComputeGravity(int stepNumber) {
  idfx::pushRegion("Gravity::ComputeGravity");
  if(havePotential) {
    if(haveUserDefPotential) {
      if(gravPotentialFunc == nullptr) {
        IDEFIX_ERROR("Gravitational potential is enabled, "
                   "but no user-defined potential has been enrolled.");
      }
      idfx::pushRegion("Gravity::user-defined:gravPotentialFunc");
      gravPotentialFunc(*data, data->t, data->x[IDIR], data->x[JDIR], data->x[KDIR], phiP);
      idfx::popRegion();
    } else {
      ResetPotential();
    }
    if(haveCentralMassPotential) {
      AddCentralMassPotential();
    }
    if(havePlanetsPotential) {
      data->planetarySystem->AddPlanetsPotential(phiP, data->t);
    }
    if(haveSelfGravityPotential) {
      // Solving Poisson for the current gas density distribution
      if(stepNumber % selfGravity.skipSelfGravity == 0) selfGravity.SolvePoisson();

      // Adding gas self-gravity contribution to global gravity potential
      selfGravity.AddSelfGravityPotential(phiP);
    }
  }
  if(haveBodyForce) {
    // When bodyforce is enabled, it has to be a userdefined function, enrolled previously
    if(bodyForceFunc == nullptr) {
      IDEFIX_ERROR("Gravitational potential is enabled, "
                   "but no user-defined potential has been enrolled.");
    }
    idfx::pushRegion("Gravity: user-defined:bodyForceFunc");
    bodyForceFunc(*data, data->t, bodyForceVector);
    idfx::popRegion();
  }

  // For debug purpose
  #ifdef DEBUG_GRAVITY
    potTotFile.open("potTot.dat",std::ios::trunc);
    this->selfGravity.WriteField(potTotFile,this->phiP);
    potTotFile.close();
  #endif

  idfx::popRegion();
}

void Gravity::EnrollPotential(GravPotentialFunc myFunc) {
  if(!this->haveUserDefPotential) {
    IDEFIX_WARNING("In order to enroll your gravitational potential, "
                 "you need to enable it first in the .ini file "
                 "with the potential entry in [Gravity].");
  }
  this->gravPotentialFunc = myFunc;
}

void Gravity::EnrollBodyForce(BodyForceFunc myFunc) {
  if(!this->haveBodyForce) {
    IDEFIX_WARNING("In order to enroll your body force, "
                 "you need to enable it first in the .ini file "
                 "with the bodyForce entry in [Gravity].");
  }
  this->bodyForceFunc = myFunc;
}

// Fill the gravitational potential with zeros
void Gravity::ResetPotential() {
  idfx::pushRegion("Gravity::ResetPotential");
  IdefixArray3D<real> phiP = this->phiP;
  idefix_for("Gravity::ResetPotential",
              0, data->np_tot[KDIR],
              0, data->np_tot[JDIR],
              0, data->np_tot[IDIR],
              KOKKOS_LAMBDA(int k, int j, int i) {
                phiP(k,j,i) = ZERO_F;
              });
  idfx::popRegion();
}

void Gravity::AddCentralMassPotential() {
  idfx::pushRegion("Gravity::AddCentralMassPotential");
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  IdefixArray1D<real> x3 = data->x[KDIR];
  IdefixArray3D<real> phiP = this->phiP;
  real mass = this->centralMass;
  real gravCst = this->gravCst;
  #if GEOMETRY == CARTESIAN
    IDEFIX_ERROR("Central mass potential is not defined in cartesian geometry");
  #endif
  idefix_for("Gravity::AddCentralMassPotential",
              0, data->np_tot[KDIR],
              0, data->np_tot[JDIR],
              0, data->np_tot[IDIR],
              KOKKOS_LAMBDA(int k, int j, int i) {
                real r;
                #if GEOMETRY == POLAR
                  r = D_EXPAND( x1(i)*x1(i),                , + x3(k)*x3(k));
                  r = sqrt(r);
                #elif GEOMETRY == CYLINDRICAL
                  r = D_EXPAND( x1(i)*x1(i), + x2(j)*x2(j)  ,              );
                  r = sqrt(r);
                #elif GEOMETRY == SPHERICAL
                  r = x1(i);
                #else
                  r = ONE_F; // Make sure this one is initialized
                #endif
                  phiP(k,j,i) += -gravCst*mass/r;
              });
  idfx::popRegion();
}
