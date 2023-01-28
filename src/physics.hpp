// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef PHYSICS_HPP_
#define PHYSICS_HPP_

#include "idefix.hpp"

// Physics type
struct Physics {
  #if HAVE_ENERGY == 1
    static constexpr bool pressure{true};
  #else
    static constexpr bool pressure{false};
  #endif
  #if MHD == YES
    static constexpr bool mhd{true};
    static constexpr int nvar{1+2*COMPONENTS  + (pressure?1:0)};
  #else
    static constexpr bool mhd{false};
    static constexpr int nvar{1+COMPONENTS  + (pressure?1:0)};
  #endif


  // prefix
  static constexpr std::string_view prefix = "Hydro";
};

#endif // PHYSICS_HPP_
