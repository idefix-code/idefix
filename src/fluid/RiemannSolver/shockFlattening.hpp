// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_RIEMANNSOLVER_SHOCKFLATTENING_HPP_
#define FLUID_RIEMANNSOLVER_SHOCKFLATTENING_HPP_

#include "fluid.hpp"
#include "../physics.hpp"
using Hydro = Fluid<Physics>;

enum class FlagShock{None, Shock};

class ShockFlattening {
 public:
  ShockFlattening(Hydro*, real);
  ShockFlattening() {}

  void FindShock();

  Hydro *hydro;
  IdefixArray3D<FlagShock> flagArray;
  bool isActive{false};
  real smoothing{0};
};
#endif // FLUID_RIEMANNSOLVER_SHOCKFLATTENING_HPP_
