// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>

#include "gravity.hpp"
#include "dataBlock.hpp"
#include "input.hpp"

void Gravity::Init(Input &input, DataBlock *datain) {
  this->data = datain;
  data->haveGravity = true;
  // Gravitational potential
  int nPotential = input.CheckEntry("Gravity","potential");
  if(nPotential >=0) {
    this->havePotential = true;
    for(int i = 0 ; i <= nPotential ; i++) {
      std::string potentialString = input.GetString("Gravity","potential",i);
      if(potentialString.compare("userdef") == 0) {
        this->haveUserDefPotential = true;
        idfx::cout << "Gravity: Enabling user-defined gravitational potential" << std::endl;
      } else if (potentialString.compare("central")) {
        this->haveCentralMassPotential = true;
        idfx::cout << "Gravity: Enabling central mass gravitational potential" << std::endl;
      } else if (potentialString.compare("selfgravity")) {
        this->haveSelfGravityPotential = true;
        idfx::cout << "Gravity: Enabling self Gravity" << std::endl;
      } else {
        IDEFIX_ERROR("Unknown type of gravitational potential in idefix.ini. ");
      }
    }
  }

  // Body Force
  if(input.CheckEntry("Gravity","bodyForce")>=0) {
    std::string potentialString = input.GetString("Gravity","bodyForce",0);
    if(potentialString.compare("userdef") == 0) {
      this->haveBodyForce = true;
      idfx::cout << "Gravity:: Enabling user-defined body force" << std::endl;
    } else {
      IDEFIX_ERROR("Unknown type of body force in idefix.ini. "
                   "Only userdef is implemented");
    }
  }

  // Allocate required arrays
  if(havePotential && !haveInitialisedPotential) {
    phiP = IdefixArray3D<real>("Gravity_PhiP",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
    haveInitialisedPotential = true;
  }
  if(haveBodyForce && !haveInitialisedBodyForce) {
    bodyForceVector = IdefixArray4D<real>("Gravity_bodyForce", COMPONENTS,
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
    haveInitialisedBodyForce = true;
  }
}

// This function compute the gravitational field, using both body force and potential
void Gravity::ComputeGravity() {
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
    // todo: implement central mass potential
    if(havePlanetsPotential) {
      IDEFIX_ERROR("Planet potential not implemented. Ask GWF.");
    }
    if(haveSelfGravityPotential) {
      IDEFIX_ERROR("Self gravity potential not implemented. Ask JM.");
    }
  }
  if(haveBodyForce) {
    // When bodyforce is enabled, it has to be a userdefined function, enrolled previously
    if(bodyForceFunc == nullptr) {
      IDEFIX_ERROR("Gravitational potential is enabled, "
                   "but no user-defined potential has been enrolled.");
    }
    idfx::pushRegion("Gravity::user-defined:bodyForceFunc");
    bodyForceFunc(*data, data->t, bodyForceVector);
    idfx::popRegion();
  }
  idfx::popRegion();
}

void Gravity::EnrollPotential(GravPotentialFunc myFunc) {
  if(!this->haveUserDefPotential) {
    IDEFIX_ERROR("In order to enroll your gravitational potential, "
                 "you need to enable it first in the .ini file "
                 "with the potential entry in [Gravity].");
  }
  this->gravPotentialFunc = myFunc;
  idfx::cout << "Gravity: User-defined gravitational potential has been enrolled" << std::endl;
}

void Gravity::EnrollBodyForce(BodyForceFunc myFunc) {
  if(!this->haveBodyForce) {
    IDEFIX_ERROR("In order to enroll your body force, "
                 "you need to enable it first in the .ini file "
                 "with the bodyForce entry in [Gravity].");
  }
  this->bodyForceFunc = myFunc;
  idfx::cout << "Gravity: User-defined body force function has been enrolled" << std::endl;
}

// Fill the gravitational potential with zeros
void Gravity::ResetPotential() {
  idfx::pushRegion("Gravity::ResetPotential");
  idefix_for("Gravity::ResetPotential",
              0, data->np_tot[KDIR],
              0, data->np_tot[JDIR],
              0, data->np_tot[IDIR],
              KOKKOS_LAMBDA(int k, int j, int i) {
                phiP(k,j,i) = ZERO_F;
              });
  idfx::popRegion();
}
