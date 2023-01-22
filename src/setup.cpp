// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "setup.hpp"
#include "idefix.hpp"

// Default setup functions. These are automatically redefined if the user provide her/his
// own implementation of the constructor, initflow and destructor

__attribute__((weak)) Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  IDEFIX_WARNING("Caution, this is the default Setup constructor and it does nothing!");
}

__attribute__((weak)) void Setup::InitFlow(DataBlock &data) {
  IDEFIX_ERROR("Please create your own setup.cpp following the documentation");
}

__attribute__((weak)) Setup::~Setup() {
}
