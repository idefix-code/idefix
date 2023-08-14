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
// The default Physics class
struct DefaultPhysics {
  static constexpr bool dust{false};
  #if HAVE_ENERGY == 1
    static constexpr bool pressure{true};
  #else
    static constexpr bool pressure{false};
  #endif
  #ifdef ISOTHERMAL
  static constexpr bool isothermal{true};
  #else
  static constexpr bool isothermal{false};
  #endif
  static constexpr bool eos = isothermal || pressure; // whether we have a eos
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

// Some Dust
struct DustPhysics {
  static constexpr bool dust{true};
  static constexpr bool pressure{false};
  static constexpr bool isothermal{false};
  static constexpr bool eos{false};

  static constexpr bool mhd{false};
  static constexpr int nvar{1+COMPONENTS};

  // prefix
  static constexpr std::string_view prefix = "Dust";
};

#endif // PHYSICS_HPP_
