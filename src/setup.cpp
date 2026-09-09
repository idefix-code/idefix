// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/setup.cpp
//
// Last modified : 03/2026
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG / CNRS - 2021 - 2026)
// - and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "setup.hpp"
#include "idefix.hpp"

// Default setup functions. These are automatically redefined if the user provide her/his
// own implementation of the constructor, initflow and destructor

__attribute__((weak)) Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  #ifndef WITH_PYTHON
    IDEFIX_WARNING("Caution, this is the default Setup constructor and it does nothing!");
  #endif
}

__attribute__((weak)) void Setup::InitFlow(DataBlock &data) {
  IDEFIX_ERROR("Please create your own setup.cpp following the documentation");
}

__attribute__((weak)) Setup::~Setup() {
}
