// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_THERMALDIFFUSION_HPP_
#define HYDRO_THERMALDIFFUSION_HPP_

#include "idefix.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "hydro_defs.hpp"


// Forward class hydro declaration
class Hydro;
class DataBlock;

class ThermalDiffusion {
 public:
  ThermalDiffusion();  // Default (empty) constructor

  void Init(Input &, Grid &, Hydro *);  // Initialisation

  void ShowConfig(); // display configuration

  void AddDiffusiveFlux(int, const real);

  // Enroll user-defined viscous diffusivity
  void EnrollThermalDiffusivity(DiffusivityFunc);

  IdefixArray4D<real> viscSrc;  // Source terms of the viscous operator
  IdefixArray3D<real> kappaArr;

  // pre-computed geometrical factors in non-cartesian geometry
  IdefixArray1D<real> one_dmu;

 private:
  Hydro *hydro; // My parent hydro object

  // type of viscosity function
  HydroModuleStatus haveThermalDiffusion{Disabled};
  DiffusivityFunc diffusivityFunc;

  // constant diffusion coefficient (when needed)
  real kappa;
};

#endif // HYDRO_THERMALDIFFUSION_HPP_
