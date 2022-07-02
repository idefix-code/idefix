// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_SHOCKFLATTENING_HPP_
#define HYDRO_SHOCKFLATTENING_HPP_

#include "hydro.hpp"

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
#endif //HYDRO_SHOCKFLATTENING_HPP_
