// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/fluid/fluid_defs.hpp
//
// Last modified : 01/2025
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2020 - 2023)
// - Clément Robert <clement.robert@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2021)
// - Victor Réville <vreville@irap.omp.eu> (IRAP - 2024 - 2025)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_FLUID_DEFS_HPP_
#define FLUID_FLUID_DEFS_HPP_

#include <vector>
#include <functional>
#include "../idefix.hpp"


// Common definitions for all of the objects dependent on hydro

// forward class declaration
class DataBlock;

template<typename Phys>
class Fluid;

#define     SMALL_PRESSURE_FIX      (1.0e-5)
#define     eps_UCT_CONTACT         (1.0e-6)


// Parabolic terms can have different status
enum HydroModuleStatus {Disabled, Constant, UserDefFunction};

// Structure to describe the status of parabolic modules
struct ParabolicModuleStatus {
  HydroModuleStatus status{Disabled};
  bool isExplicit{false};
  bool isRKL{false};
};


using GravPotentialFunc = std::function<void(DataBlock &, const real, IdefixArray1D<real>&,
                                     IdefixArray1D<real>&, IdefixArray1D<real>&,
                                     IdefixArray3D<real> &)>;

using BodyForceFunc = std::function<void(DataBlock &, const real, IdefixArray4D<real>&)>;

using IsoSoundSpeedFunc = std::function<void(DataBlock &, const real, IdefixArray3D<real> &)>;

template<typename Phys>
using SrcTermFunc = std::function<void(Fluid<Phys> *, const real, const real)>;


using EmfBoundaryFunc = std::function<void(DataBlock &, const real)>;
using DiffusivityFunc = std::function<void(DataBlock &, const real, IdefixArray3D<real> &)>;

using BragDiffusivityFunc = std::function<void(DataBlock &, const real,
                       std::vector<IdefixArray3D<real>> &)>;

// Deprecated signatures (kept as aliases for backward compatibility)
using SrcTermFuncOld = std::function<void(DataBlock &, const real, const real)>;


#endif //FLUID_FLUID_DEFS_HPP_
