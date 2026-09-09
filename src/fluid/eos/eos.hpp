// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/fluid/eos/eos.hpp
//
// Last modified : 10/2023
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG / CNRS - 2022 - 2023)
// - and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_EOS_EOS_HPP_
#define FLUID_EOS_EOS_HPP_

// This is a wrapper which decides which eos .hpp file should be chosen according to
// Idefix configuration

#ifndef EOS_FILE
  #ifdef ISOTHERMAL
    #include "eos_isothermal.hpp"
  #else
    #include "eos_adiabatic.hpp"
  #endif
#else
  #include EOS_FILE
#endif

#endif // FLUID_EOS_EOS_HPP_
