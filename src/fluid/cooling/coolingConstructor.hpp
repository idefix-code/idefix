// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_COOLING_COOLINGCONSTRUCTOR_HPP_
#define FLUID_COOLING_COOLINGCONSTRUCTOR_HPP_

#include <string>
#include "cooling.hpp"
#include "fluid.hpp"

template<typename Phys>
RadCooling::RadCooling(Input &input, Grid &grid, Fluid<Phys> *hydroin):
                       Vc(hydroin->Vc),
                       Vs(hydroin->Vs) {
  idfx::pushRegion("RadCooling::RadCooling");
  // Save the parent hydro object
  this->data = hydroin->data;
  this->eos = *(hydroin->eos.get());

  if(input.CheckEntry("Hydro","Cooling")>=0) {
    idfx::cout << input.Get<std::string>("Hydro","Cooling",1) << std::endl;
    if(input.Get<std::string>("Hydro","Cooling",0).compare("Tabulated") == 0) {
      this->cooling_type = "tabulated";
      this->cooltable_filename = input.Get<std::string>("Hydro","Cooling",1);
      if(input.Get<std::string>("Hydro","Cooling",2).compare("Townsend") == 0) {
        this->cooling_integration = "townsend";
        this->rad_cooling_int_func = [this](real dt){ TownsendIntegration(dt); };
      } else {
        IDEFIX_ERROR("Unknown radiative cooling integration in idefix.ini. "
                     "Can only be Townsend at the moment.");
      }
      if(input.Get<std::string>("Hydro","Cooling",3).compare("TcoolFloor") == 0) {
        this->TcoolFloor = input.Get<real>("Hydro","Cooling",4);
      } else {
        IDEFIX_WARNING("No cooling floor set "
                       "Defaulting to 1.0e+04 K");
        this->TcoolFloor = 1.0e+04;
      }
    /*
    else if(input.Get<std::string>("Hydro","Cooling",1).compare("neq") == 0) {
    } else if (input.Get<std::string>("Hydro","Cooling",1).compare("nolimiter") == 0) {
    */
    } else {
      IDEFIX_ERROR("Unknown radiative cooling in idefix.ini. "
                   "Can only be Tabulated at the moment.");
    }
  } else {
    IDEFIX_ERROR("I cannot create a RadCooling object without Cooling defined"
                   "in the .ini file");
  }

  InitArrays();

  idfx::popRegion();
}

#endif // FLUID_COOLING_COOLINGCONSTRUCTOR_HPP_
