// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/fluid/RiemannSolver/shockFlag.hpp
//
// Last modified : 10/2026
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2023)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_RIEMANNSOLVER_SHOCKFLAG_HPP_
#define FLUID_RIEMANNSOLVER_SHOCKFLAG_HPP_

#include <functional>

#include "arrays.hpp"

// Forward declaration (used in the enrollment function signature)
class DataBlock;

// Flag attached to each cell: whether shock flattening should be applied there
enum class FlagShock{None, Shock};

// User-defined function used to flag additional cells for shock flattening.
// It is called on the host and receives the datablock, the current time and
// the (writable) shock flag array.
using UserShockFunc = std::function<void(DataBlock &, const real, IdefixArray3D<FlagShock>&)>;

#endif // FLUID_RIEMANNSOLVER_SHOCKFLAG_HPP_
