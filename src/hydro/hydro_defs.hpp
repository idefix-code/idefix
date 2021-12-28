// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_HYDRO_DEFS_HPP_
#define HYDRO_HYDRO_DEFS_HPP_

#include "../idefix.hpp"

// Common definitions for all of the objects dependent on hydro

// forward class declaration
class DataBlock;

#define     SMALL_PRESSURE_FIX      (1.0e-5)
#define     eps_UCT_CONTACT         (1.0e-6)

// Riemann Solver type
#if MHD == YES
enum Solver {TVDLF=1, HLL, HLLD, ROE};
#else
enum Solver {TVDLF=1, HLL, HLLC, ROE};
#endif

// Parabolic terms can have different status
enum HydroModuleStatus {Disabled, Constant, UserDefFunction };

// Structure to describe the status of parabolic modules
struct ParabolicModuleStatus {
  HydroModuleStatus status{Disabled};
  bool isExplicit{false};
  bool isRKL{false};
};


using UserDefBoundaryFunc = void (*) (DataBlock &, int dir, BoundarySide side,
                                      const real t);
using GravPotentialFunc = void (*) (DataBlock &, const real t, IdefixArray1D<real>&,
                                    IdefixArray1D<real>&, IdefixArray1D<real>&,
                                    IdefixArray3D<real> &);

using BodyForceFunc = void (*) (DataBlock &, const real t, IdefixArray4D<real>&);

using SrcTermFunc = void (*) (DataBlock &, const real t, const real dt);
using InternalBoundaryFunc = void (*) (DataBlock &, const real t);
using EmfBoundaryFunc = void (*) (DataBlock &, const real t);
using DiffusivityFunc = void (*) (DataBlock &, const real t, IdefixArray3D<real> &);
using IsoSoundSpeedFunc = void (*) (DataBlock &, const real t, IdefixArray3D<real> &);

#endif //HYDRO_HYDRO_DEFS_HPP_
