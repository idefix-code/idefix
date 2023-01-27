// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_THERMALDIFFUSION_HPP_
#define FLUID_THERMALDIFFUSION_HPP_

#include <string>

#include "idefix.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "fluid_defs.hpp"


// Forward class hydro declaration
template <typename Phys> class Fluid;

class DataBlock;

class ThermalDiffusion {
 public:
  template <typename Phys>
  ThermalDiffusion(Input &, Grid &, Fluid<Phys> *);  // Initialisation

  void ShowConfig(); // display configuration

  void AddDiffusiveFlux(int, const real, const IdefixArray4D<real> &);

  // Enroll user-defined viscous diffusivity
  void EnrollThermalDiffusivity(DiffusivityFunc);

  IdefixArray4D<real> viscSrc;  // Source terms of the viscous operator
  IdefixArray3D<real> kappaArr;

  // pre-computed geometrical factors in non-cartesian geometry
  IdefixArray1D<real> one_dmu;

 private:
  DataBlock *data;

  // status of the module
  ParabolicModuleStatus &status;

  DiffusivityFunc diffusivityFunc;

  // helper array
  IdefixArray4D<real> &Vc;
  IdefixArray3D<real> &dMax;

  // constant diffusion coefficient (when needed)
  real kappa;

  // adiabatic exponent (required to get the heat capacity)
  // TODO(glesur): generalize to any equation of state!
  real gamma;
};

#include "fluid.hpp"
#include "dataBlock.hpp"

template <typename Phys>
ThermalDiffusion::ThermalDiffusion(Input &input, Grid &grid, Fluid<Phys> *hydroin):
                            Vc{hydroin->Vc},
                            dMax{hydroin->dMax},
                            gamma{hydroin->GetGamma()},
                            data{hydroin->data},
                            status{hydroin->thermalDiffusionStatus} {
  idfx::pushRegion("ThermalDiffusion::ThermalDiffusion");

  if(input.CheckEntry("Hydro","TDiffusion")>=0) {
    if(input.Get<std::string>("Hydro","TDiffusion",1).compare("constant") == 0) {
        this->kappa = input.Get<real>("Hydro","TDiffusion",2);
        this->status.status = Constant;
      } else if(input.Get<std::string>("Hydro","TDiffusion",1).compare("userdef") == 0) {
        this->status.status = UserDefFunction;
        this->kappaArr = IdefixArray3D<real>("ThermalDiffusionKappaArray",data->np_tot[KDIR],
                                                                 data->np_tot[JDIR],
                                                                 data->np_tot[IDIR]);

      } else {
        IDEFIX_ERROR("Unknown thermal diffusion definition in idefix.ini. "
                     "Can only be constant or userdef.");
      }
  } else {
    IDEFIX_ERROR("I cannot create a ThermalDiffusion object without TDiffusion defined"
                   "in the .ini file");
  }

  #ifdef ISOTHERMAL
    IDEFIX_ERROR("Thermal diffusion is not compatible with the ISOTHERMAL approximation");
  #endif

  idfx::popRegion();
}

#endif // FLUID_THERMALDIFFUSION_HPP_
