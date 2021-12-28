// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "setup.hpp"
#include "idefix.hpp"


Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  IDEFIX_WARNING("Caution, this is the default Setup constructor and it does nothing!");
}

void Setup::InitFlow(DataBlock &data) {
  IDEFIX_ERROR("Please create your own setup.cpp following the documentation");
}
