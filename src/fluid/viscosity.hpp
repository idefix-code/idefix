// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_VISCOSITY_HPP_
#define FLUID_VISCOSITY_HPP_

#include <string>

#include "idefix.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "fluid_defs.hpp"



// Forward class hydro declaration
template <typename Phys> class Fluid;
class DataBlock;

using ViscousDiffusivityFunc = void (*) (DataBlock &, const real t,
                                         IdefixArray3D<real> &, IdefixArray3D<real> &);

class Viscosity {
 public:
  template <typename Phys>
  Viscosity(Input &, Grid &, Fluid<Phys> *);
  void ShowConfig();                    // print configuration
  void AddViscousFlux(int, const real, const IdefixArray4D<real> &);

  // Enroll user-defined viscous diffusivity
  void EnrollViscousDiffusivity(ViscousDiffusivityFunc);

  // Function for internal use (but public to allow for Cuda lambda capture)
  void InitArrays();

  IdefixArray4D<real> viscSrc;  // Source terms of the viscous operator
  IdefixArray3D<real> eta1Arr;
  IdefixArray3D<real> eta2Arr;

  // pre-computed geometrical factors in non-cartesian geometry
  IdefixArray1D<real> one_dmu;

 private:
  DataBlock* data;

  // Viscosity status
  ParabolicModuleStatus &status;

  ViscousDiffusivityFunc viscousDiffusivityFunc;

  IdefixArray4D<real> &Vc;
  IdefixArray3D<real> &dMax;

  // constant diffusion coefficient (when needed)
  real eta1, eta2;

  // Shearing box
  real sbS;
};

#include "fluid.hpp"

template<typename Phys>
Viscosity::Viscosity(Input &input, Grid &grid, Fluid<Phys> *hydroin):
                      Vc{hydroin->Vc},
                      dMax{hydroin->dMax},
                      status{hydroin->viscosityStatus} {
  idfx::pushRegion("Viscosity::Viscosity");
  // Save the parent hydro object
  this->data = hydroin->data;
  this->sbS = hydroin->sbS;

  if(status.status == Constant) {
    this->eta1 = input.Get<real>(std::string(Phys::prefix),"viscosity",2);
    // second viscosity?
    this->eta2 = input.GetOrSet<real>(std::string(Phys::prefix),"viscosity",3, 0.0);
  } else if(status.status == UserDefFunction) {
    this->eta1Arr = IdefixArray3D<real>("ViscosityEta1Array",data->np_tot[KDIR],
                                                              data->np_tot[JDIR],
                                                              data->np_tot[IDIR]);
    this->eta2Arr = IdefixArray3D<real>("ViscosityEta1Array",data->np_tot[KDIR],
                                                              data->np_tot[JDIR],
                                                              data->np_tot[IDIR]);
  } else {
        IDEFIX_ERROR("Unknown viscosity definition in idefix.ini. "
                     "Can only be constant or userdef.");
  }

  InitArrays();

  idfx::popRegion();
}
#endif // FLUID_VISCOSITY_HPP_
