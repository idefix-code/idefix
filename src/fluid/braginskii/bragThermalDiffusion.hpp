// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_BRAGTHERMALDIFFUSION_HPP_
#define FLUID_BRAGTHERMALDIFFUSION_HPP_

#include <string>

#include "idefix.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "fluid_defs.hpp"

real minmodTh(const real, const real); 
real vanLeerTh(const real, const real);
real monotonizedCentralTh(const real, const real); 

// Forward class hydro declaration
template <typename Phys> class Fluid;

class DataBlock;

class BragThermalDiffusion {
 public:
  template <typename Phys>
  BragThermalDiffusion(Input &, Grid &, Fluid<Phys> *);  // Initialisation

  void ShowConfig(); // display configuration

  void AddDiffusiveFlux(int, const real, const IdefixArray4D<real> &);

  // Enroll user-defined thermal conductivity
  void EnrollBragThermalDiffusivity(BragDiffusivityFunc);

  IdefixArray3D<real> heatSrc;  // Source terms of the thermal operator
  IdefixArray3D<real> knorArr;
  IdefixArray3D<real> kparArr;

  // pre-computed geometrical factors in non-cartesian geometry
  IdefixArray1D<real> one_dmu;

 private:
  DataBlock *data;

  // status of the module
  ParabolicModuleStatus &status;

  BragDiffusivityFunc diffusivityFunc;

  bool haveSlopeLimiter;

  // helper array
  IdefixArray4D<real> &Vc;
  IdefixArray4D<real> &Vs;
  IdefixArray3D<real> &dMax;

  // constant diffusion coefficient (when needed)
  real knor, kpar;

  // adiabatic exponent (required to get the heat capacity)
  // TODO(glesur): generalize to any equation of state!
  real gamma;

//  // type of thermal conductivity function
//  bool haveBraginskiiConductivity;
//  bool haveSlopeLimiter;
//  HydroModuleStatus haveThConductivity;
//  ThermalConductivityFunc thermalConductivityFunc;
//  BragThermalConductivityFunc bragThermalConductivityFunc;
  SlopeLimiterFunc slopeLimiter;
};

#include "fluid.hpp"
#include "dataBlock.hpp"

template <typename Phys>
BragThermalDiffusion::BragThermalDiffusion(Input &input, Grid &grid, Fluid<Phys> *hydroin):
                            Vc{hydroin->Vc},
                            Vs{hydroin->Vs},
                            dMax{hydroin->dMax},
                            gamma{hydroin->GetGamma()},
                            data{hydroin->data},
                            status{hydroin->bragThermalDiffusionStatus} {
  idfx::pushRegion("BragThermalDiffusion::BragThermalDiffusion");

  if(input.CheckEntry("Hydro","bragTDiffusion")>=0) {
    if(input.Get<std::string>("Hydro","bragTDiffusion",1).compare("constant") == 0) {
      this->haveSlopeLimiter = false;
      this->kpar = input.Get<real>("Hydro","bragTDiffusion",2);
      this->knor = input.GetOrSet<real>("Hydro","bragTDiffusion",3,0.);
      this->status.status = Constant;
    } else if(input.Get<std::string>("Hydro","bragTDiffusion",1).compare("userdef") == 0) {
      this->haveSlopeLimiter = false;
      this->status.status = UserDefFunction;
      this->kparArr = IdefixArray3D<real>("BragThermalDiffusionKparArray",data->np_tot[KDIR],
                                                               data->np_tot[JDIR],
                                                               data->np_tot[IDIR]);
      this->knorArr = IdefixArray3D<real>("BragThermalDiffusionKnorArray",data->np_tot[KDIR],
                                                               data->np_tot[JDIR],
                                                               data->np_tot[IDIR]);
    } else if (input.Get<std::string>("Hydro","bragTDiffusion",1).compare("limiter") == 0) {
      this->haveSlopeLimiter = true;
      if(input.Get<std::string>("Hydro","bragTDiffusion",2).compare("minmod") == 0) {
        slopeLimiter = minmodTh;
      } else if(input.Get<std::string>("Hydro","bragTDiffusion",2).compare("vanleer") == 0) {
        slopeLimiter = vanLeerTh;
      } else if(input.Get<std::string>("Hydro","bragTDiffusion",2).compare("mc") == 0) {
        slopeLimiter = monotonizedCentralTh;
      } else {
        IDEFIX_ERROR("Unknown braginskii thermal diffusion limiter in idefix.ini. "
                     "Can only be minmod, vanleer or mc.");
      }
      if(input.Get<std::string>("Hydro","bragTDiffusion",3).compare("constant") == 0) {
          this->kpar = input.Get<real>("Hydro","bragTDiffusion",4);
          this->knor = input.Get<real>("Hydro","bragTDiffusion",4);
          this->status.status = Constant;
      } else if(input.Get<std::string>("Hydro","bragTDiffusion",3).compare("userdef") == 0) {
          this->status.status = UserDefFunction;
          this->kparArr = IdefixArray3D<real>("BragThermalDiffusionKparArray",data->np_tot[KDIR],
                                                                   data->np_tot[JDIR],
                                                                   data->np_tot[IDIR]);
          this->knorArr = IdefixArray3D<real>("BragThermalDiffusionKnorArray",data->np_tot[KDIR],
                                                                   data->np_tot[JDIR],
                                                                   data->np_tot[IDIR]);
      } else {
        IDEFIX_ERROR("Unknown braginskii thermal diffusion definition in idefix.ini. "
                     "Can only be constant or userdef.");
      }
    } else {
      IDEFIX_ERROR("Unknown braginskii thermal diffusion definition in idefix.ini. "
                   "Can only be constant or userdef.");
    }
  } else {
    IDEFIX_ERROR("I cannot create a BragThermalDiffusion object without bragTDiffusion defined"
                   "in the .ini file");
  }

  #ifndef MHD
    IDEFIX_ERROR("Braginskii Thermal diffusion requires MHD");
  #endif
  #ifdef ISOTHERMAL
    IDEFIX_ERROR("Braginskii Thermal diffusion is not compatible with the ISOTHERMAL approximation");
  #endif

  idfx::popRegion();
}

//using ThermalConductivityFunc = void (*) (DataBlock &, const real t, IdefixArray3D<real> &);
//using BragThermalConductivityFunc = void (*) (DataBlock &, const real t,
//                                        IdefixArray3D<real> &, IdefixArray3D<real> &);
//                                        // normal and parralel conductivities
//using SlopeLimiterFunc = real (*) (const real, const real);

#endif // FLUID_BRAGTHERMALDIFFUSION_HPP_
