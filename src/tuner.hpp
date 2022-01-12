// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef TUNER_HPP_
#define TUNER_HPP_

#include "idefix.hpp"
#include "dataBlock.hpp"
#include "input.hpp"
#include "setup.hpp"

namespace Tuner {
void tuneLoops(DataBlock &, Setup &, Input &);
}

#endif // TUNER_HPP_
