// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_BRAGINSKII_BRAGVISCOSITY_HPP_
#define FLUID_BRAGINSKII_BRAGVISCOSITY_HPP_

#include <string>

#include "idefix.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "fluid_defs.hpp"

real minmodV(const real, const real);
real vanLeerV(const real, const real);
real monotonizedCentralV(const real, const real);

// Forward class hydro declaration
template <typename Phys> class Fluid;
class DataBlock;

class BragViscosity {
 public:
  template <typename Phys>
  BragViscosity(Input &, Grid &, Fluid<Phys> *);
  void ShowConfig();                    // print configuration
  void AddBragViscousFlux(int, const real, const IdefixArray4D<real> &);

  // Enroll user-defined viscous diffusivity
  void EnrollBragViscousDiffusivity(DiffusivityFunc);

  // Function for internal use (but public to allow for Cuda lambda capture)
  void InitArrays();

  IdefixArray4D<real> bragViscSrc;  // Source terms of the viscous operator
  IdefixArray3D<real> etaBragArr;


  // pre-computed geometrical factors in non-cartesian geometry
  IdefixArray1D<real> one_dmu;

 private:
  DataBlock* data;

  // Viscosity status
  ParabolicModuleStatus &status;

  DiffusivityFunc bragViscousDiffusivityFunc;

  bool haveSlopeLimiter;

  IdefixArray4D<real> &Vc;
  IdefixArray4D<real> &Vs;
  IdefixArray3D<real> &dMax;

  // constant diffusion coefficient (when needed)
  real etaBrag;

  // Shearing box
  real sbS;
//  // type of viscosity function
//  bool haveBraginskiiViscosity;
//  bool haveSlopeLimiter;
//  HydroModuleStatus haveViscosity{Disabled};
//  ViscousDiffusivityFunc viscousDiffusivityFunc;
//  BragViscousDiffusivityFunc bragViscousDiffusivityFunc;
  SlopeLimiterFunc slopeLimiter;
};

#include "fluid.hpp"

template<typename Phys>
BragViscosity::BragViscosity(Input &input, Grid &grid, Fluid<Phys> *hydroin):
                      Vc{hydroin->Vc},
                      Vs{hydroin->Vs},
                      dMax{hydroin->dMax},
                      status{hydroin->bragViscosityStatus} {
  idfx::pushRegion("BragViscosity::BragViscosity");
  // Save the parent hydro object
  this->data = hydroin->data;
  this->sbS = hydroin->sbS;

  if(input.CheckEntry("Hydro","bragViscosity")>=0) {
    if(input.Get<std::string>("Hydro","bragViscosity",1).compare("constant") == 0) {
        this->etaBrag = input.Get<real>("Hydro","bragViscosity",2);
        this->status.status = Constant;
    } else if(input.Get<std::string>("Hydro","bragViscosity",1).compare("userdef") == 0) {
        this->status.status = UserDefFunction;
        this->etaBragArr = IdefixArray3D<real>("BragViscosityEtaArray",data->np_tot[KDIR],
                                                                 data->np_tot[JDIR],
                                                                 data->np_tot[IDIR]);
    } else if (input.Get<std::string>("Hydro","bragViscosity",1).compare("limiter") == 0) {
      this->haveSlopeLimiter = true;
      if(input.Get<std::string>("Hydro","bragViscosity",2).compare("minmod") == 0) {
        slopeLimiter = minmodV;
      } else if(input.Get<std::string>("Hydro","bragViscosity",2).compare("vanleer") == 0) {
        slopeLimiter = vanLeerV;
      } else if(input.Get<std::string>("Hydro","bragViscosity",2).compare("mc") == 0) {
        slopeLimiter = monotonizedCentralV;
      } else {
        IDEFIX_ERROR("Unknown braginskii viscosity limiter in idefix.ini. "
                     "Can only be minmod, vanleer or mc.");
      }
      if(input.Get<std::string>("Hydro","bragViscosity",3).compare("constant") == 0) {
          this->etaBrag = input.Get<real>("Hydro","bragViscosity",4);
          this->status.status = Constant;
      } else if(input.Get<std::string>("Hydro","bragViscosity",3).compare("userdef") == 0) {
          this->status.status = UserDefFunction;
          this->etaBragArr = IdefixArray3D<real>("BragViscosityEtaArray",data->np_tot[KDIR],
                                                                   data->np_tot[JDIR],
                                                                   data->np_tot[IDIR]);
      } else {
        IDEFIX_ERROR("Unknown braginskii viscosity definition in idefix.ini. "
                     "Can only be constant or userdef.");
      }
    } else {
      IDEFIX_ERROR("Unknown braginskii viscosity definition in idefix.ini. "
                   "Can only be constant or userdef.");
    }
  } else {
    IDEFIX_ERROR("I cannot create a BragViscosity object without viscosity defined"
                   "in the .ini file");
  }

  InitArrays();

  idfx::popRegion();
}
#endif // FLUID_BRAGINSKII_BRAGVISCOSITY_HPP_
