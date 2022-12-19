// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef PHYSICS_HPP_
#define PHYSICS_HPP_

#include "idefix.hpp"

// Physics type
struct Physics {
  static constexpr bool mhd{true};
  #if HAVE_ENERGY == 1
    static constexpr bool pressure{true};
  #else
    static constexpr bool pressure{false};
  #endif
  static constexpr int nvar{1+2*COMPONENTS  + (pressure?1:0)};
  // prefix 
  static constexpr std::string_view prefix = "Hydro";
};

#endif // PHYSICS_HPP_