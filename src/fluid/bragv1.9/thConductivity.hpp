// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_THCONDUCTIVITY_HPP_
#define HYDRO_THCONDUCTIVITY_HPP_

#include "idefix.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "hydro_defs.hpp"


// Forward class hydro declaration
class Hydro;
class DataBlock;

using ThermalConductivityFunc = void (*) (DataBlock &, const real t, IdefixArray3D<real> &);
using BragThermalConductivityFunc = void (*) (DataBlock &, const real t,
                                        IdefixArray3D<real> &, IdefixArray3D<real> &);
                                        // normal and parralel conductivities
using SlopeLimiterFunc = real (*) (const real, const real);

real minmodTh(const real, const real); 
real vanLeerTh(const real, const real);
real monotonizedCentralTh(const real, const real); 

class ThConductivity {
 public:
  ThConductivity();  // Default (empty) constructor

  void Init(Input &, Grid &, Hydro *);  // Initialisation

  void AddThermalFlux(int, const real);

  void ShowConfig();

  // Enroll user-defined thermal conductivity
  void EnrollThermalConductivity(ThermalConductivityFunc);
  void EnrollBragThermalConductivity(BragThermalConductivityFunc);

  IdefixArray3D<real> heatSrc;  // Source terms of the thermal operator

  IdefixArray3D<real> kappaArr;
  IdefixArray3D<real> knorArr;
  IdefixArray3D<real> kparArr;

  // pre-computed geometrical factors in non-cartesian geometry
  IdefixArray1D<real> one_dmu;

 private:
  Hydro *hydro; // My parent hydro object

  // type of thermal conductivity function
  bool haveBraginskiiConductivity;
  bool haveSlopeLimiter;
  HydroModuleStatus haveThConductivity;
  ThermalConductivityFunc thermalConductivityFunc;
  BragThermalConductivityFunc bragThermalConductivityFunc;
  SlopeLimiterFunc slopeLimiter;

  // constant diffusion coefficient (when needed)
  real kappa, knor, kpar;
};

#endif // HYDRO_THCONDUCTIVITY_HPP_
