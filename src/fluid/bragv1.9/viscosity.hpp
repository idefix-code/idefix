// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_VISCOSITY_HPP_
#define HYDRO_VISCOSITY_HPP_

#include "idefix.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "hydro_defs.hpp"


// Forward class hydro declaration
class Hydro;
class DataBlock;

using ViscousDiffusivityFunc = void (*) (DataBlock &, const real t,
                                         IdefixArray3D<real> &, IdefixArray3D<real> &);
using BragViscousDiffusivityFunc = void (*) (DataBlock &, const real t, IdefixArray3D<real> &);
using SlopeLimiterFunc = real (*) (const real, const real);

class Viscosity {
 public:
  Viscosity();  // Default (empty) constructor

  void Init(Input &, Grid &, Hydro *);  // Initialisation
  void ShowConfig();                    // print configuration

  void AddViscousFlux(int, const real);

  // Enroll user-defined viscous diffusivity
  void EnrollViscousDiffusivity(ViscousDiffusivityFunc);
  void EnrollBragViscousDiffusivity(BragViscousDiffusivityFunc);

  IdefixArray4D<real> viscSrc;  // Source terms of the viscous operator
  IdefixArray3D<real> eta1Arr;
  IdefixArray3D<real> eta2Arr;
  IdefixArray3D<real> etaBragArr;

  // pre-computed geometrical factors in non-cartesian geometry
  IdefixArray1D<real> one_dmu;

 private:
  Hydro *hydro; // My parent hydro object

  // type of viscosity function
  bool haveBraginskiiViscosity;
  bool haveSlopeLimiter;
  HydroModuleStatus haveViscosity{Disabled};
  ViscousDiffusivityFunc viscousDiffusivityFunc;
  BragViscousDiffusivityFunc bragViscousDiffusivityFunc;
  SlopeLimiterFunc slopeLimiter;

  // constant diffusion coefficient (when needed)
  real eta1, eta2, etaBrag;
};

#endif // HYDRO_VISCOSITY_HPP_
