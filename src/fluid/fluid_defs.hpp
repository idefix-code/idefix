// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_FLUID_DEFS_HPP_
#define FLUID_FLUID_DEFS_HPP_

#include "../idefix.hpp"


// Common definitions for all of the objects dependent on hydro

// forward class declaration
class DataBlock;

template<typename Phys>
class Fluid;

#define     SMALL_PRESSURE_FIX      (1.0e-5)
#define     eps_UCT_CONTACT         (1.0e-6)



// Parabolic terms can have different status
enum HydroModuleStatus {Disabled, Constant, UserDefFunction };

// Structure to describe the status of parabolic modules
struct ParabolicModuleStatus {
  HydroModuleStatus status{Disabled};
  bool isExplicit{false};
  bool isRKL{false};
};


using GravPotentialFunc = void (*) (DataBlock &, const real t, IdefixArray1D<real>&,
                                    IdefixArray1D<real>&, IdefixArray1D<real>&,
                                    IdefixArray3D<real> &);

using BodyForceFunc = void (*) (DataBlock &, const real t, IdefixArray4D<real>&);

using IsoSoundSpeedFunc = void (*) (DataBlock &, const real t, IdefixArray3D<real> &);

template<typename Phys>
using SrcTermFunc = void (*) (Fluid<Phys> *, const real t, const real dt);


using EmfBoundaryFunc = void (*) (DataBlock &, const real t);
using DiffusivityFunc = void (*) (DataBlock &, const real t, IdefixArray3D<real> &);
using BragDiffusivityFunc = void (*) (DataBlock &, const real t,
                                      IdefixArray3D<real> &, IdefixArray3D<real> &);

// Deprecated signatures
using SrcTermFuncOld = void (*) (DataBlock &, const real t, const real dt);


#endif //FLUID_FLUID_DEFS_HPP_
