// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_TOWNSEND_RADIATIVE_HPP_
#define FLUID_TOWNSEND_RADIATIVE_HPP_

#include <string>

#include "idefix.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "fluid_defs.hpp"
#include "eos.hpp"
#include "lookupTable.hpp"

// Forward class hydro declaration
template <typename Phys> class Fluid;
class DataBlock;
class EquationOfState;

class RadCooling {
 public:
  template <typename Phys>
  RadCooling(Input &, Grid &, Fluid<Phys> *);
  void ShowConfig();                      // print configuration

  // Function for internal use (but public to allow for Cuda lambda capture)
  void InitArrays();
  void InitUnits();

  LookupTable<1> cooltable_data;
  IdefixArray1D<real> temperature_data;
  IdefixArray1D<real> Lambda_cool_data;

  real rho_unit, mass_unit, vel_unit, time_unit, len_unit;
  real TcoolFloor;

  IdefixArray3D<real> delta_eng;  // Source terms of the cooling operator
  std::string cooling_type;
  std::string cooling_integration;
  std::string cooltable_filename;

  std::function<void(real)> rad_cooling_int_func;
  void TownsendIntegration(real);
  void CalculateCoolingSource(real);            // add the source term to fluid equation

 private:
  DataBlock* data;
  EquationOfState eos;

  IdefixArray4D<real> &Vc;
  IdefixArray4D<real> &Vs;
};

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

  InitUnits();
  InitArrays();

  idfx::popRegion();
}

#endif // FLUID_TOWNSEND_RADIATIVE_HPP_
