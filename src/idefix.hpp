// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef IDEFIX_HPP_
#define IDEFIX_HPP_
#include <fstream>
#include <iostream>
#include <cstdio>
#include <Kokkos_Core.hpp>
// #include <Kokkos_DualView.hpp> // do we still need this?
#ifdef WITH_MPI
#include <mpi.h>
#endif

using Device = Kokkos::DefaultExecutionSpace;
using Layout = Kokkos::LayoutRight;

/// Type of loops we admit in idefix (see loop.hpp for details)
enum class LoopPattern { SIMDFOR, RANGE, MDRANGE, TPX, TPTTRTVR, UNDEFINED };

#define     YES     255
#define     NO      0

#define     SMALL_NUMBER      (1e-10)

#ifndef MHD
#warning MHD flag should be set to yes or no explicitly. I will assume MHD is enabled.
#define MHD  YES
#endif

/* ---- Geometry Labels ( > 0) ----  */

#define CARTESIAN    1
#define CYLINDRICAL  2
#define POLAR        3
#define SPHERICAL    4

#define     IDIR    0
#define     JDIR    1
#define     KDIR    2


// Basic configuration
#ifndef DEFINITIONS_FILE
  #include "definitions.hpp"
#else
  #include DEFINITIONS_FILE
#endif
#include "real_types.hpp"

#ifdef EMF_AVERAGE
#error EMF_AVERAGE is deprecated. Use hydro/emf in the input file to set the emf averaging scheme
#endif
#if GEOMETRY == CYLINDRICAL && DIMENSIONS == 3
  #error CYLINDRICAL should only be used with DIMENSIONS <= 2. Use POLAR for 3D problems.
#endif


// Check whether we're isothermal. If we're not, then we need to solve an energy equation
#ifndef HAVE_ENERGY
  #ifdef ISOTHERMAL
    #define HAVE_ENERGY   0
  #else
    #define HAVE_ENERGY   1
  #endif
#endif

// Shortcuts for fields used in the code

#define  RHO 0
#define  MX1 1
#define  MX2 (COMPONENTS >= 2 ? 2: 255)
#define  MX3 (COMPONENTS == 3 ? 3: 254)
#if MHD == YES
#define  BX1 (COMPONENTS + 1 + HAVE_ENERGY)
#define  BX2 (COMPONENTS >= 2 ? (BX1+1): 252)
#define  BX3 (COMPONENTS == 3 ? (BX1+2): 251)

#else
#define  BX1 253
#define  BX2 252
#define  BX3 251
#endif

#if HAVE_ENERGY == 1
  #define ENG  (COMPONENTS + 1)
#else
  #define ENG   250
#endif
#define PRS  ENG

#define VX1   MX1
#define VX2   MX2
#define VX3   MX3

#if MHD == YES
  #define NFLX (1 + 2*COMPONENTS + HAVE_ENERGY)
  #define TRG   NFLX                  // Gas Tracer index
  #define TRD  (NFLX-COMPONENTS)      // Dust Tracer index
#else
  #define NFLX (1 + COMPONENTS + HAVE_ENERGY)
  #define TRG  NFLX
  #define TRD  NFLX
#endif

// Face-centered variables
#define BX1s  0
#define BX2s  1
#define BX3s  2

// Edge-centered variables
#ifdef EVOLVE_VECTOR_POTENTIAL
  #if DIMENSIONS < 3
    #define AX1e   250
    #define AX2e   251
    #define AX3e   0
  #else
    #define AX1e   0
    #define AX2e   1
    #define AX3e   2
  #endif
#endif

// User-Friendly variables in non-cartesian geometry
#if GEOMETRY == CYLINDRICAL
  #if COMPONENTS >= 1
    #define iVR    VX1
    #define iMR    MX1
    #define iBR    BX1
  #endif

  #if COMPONENTS >= 2
    #define iVZ    VX2
    #define iMZ    MX2
    #define iBZ    BX2
  #endif

  #if COMPONENTS >= 3
    #define iVPHI  VX3
    #define iMPHI  MX3
    #define iBPHI  BX3
  #endif
#endif

#if GEOMETRY == POLAR
  #if COMPONENTS >= 1
    #define iVR    VX1
    #define iMR    MX1
    #define iBR    BX1
  #endif

  #if COMPONENTS >= 2
    #define iVPHI  VX2
    #define iMPHI  MX2
    #define iBPHI  BX2
  #endif

  #if COMPONENTS == 3
    #define iVZ    VX3
    #define iMZ    MX3
    #define iBZ    BX3
  #endif
#endif

#if GEOMETRY == SPHERICAL
  #if COMPONENTS >= 1
    #define iVR    VX1
    #define iMR    MX1
    #define iBR    BX1
  #endif

  #if COMPONENTS >= 2
    #define iVTH   VX2
    #define iMTH   MX2
    #define iBTH   BX2
  #endif

  #if COMPONENTS == 3
    #define iVPHI  VX3
    #define iMPHI  MX3
    #define iBPHI  BX3
  #endif
#endif

// Some macro definitions


#if DIMENSIONS == 1
    #define     IOFFSET     1
    #define     JOFFSET     0
    #define     KOFFSET     0
#endif
#if DIMENSIONS == 2
    #define     IOFFSET     1
    #define     JOFFSET     1
    #define     KOFFSET     0
#endif
#if DIMENSIONS == 3
    #define     IOFFSET     1
    #define     JOFFSET     1
    #define     KOFFSET     1
#endif


#define NVAR    (NFLX)


// File handler depends on the type of I/O we use
#ifdef WITH_MPI
using IdfxFileHandler = MPI_File;
#else
using IdfxFileHandler = FILE*;
#endif

// Types of boundary which can be treated
enum BoundaryType { internal, periodic, reflective, outflow, shearingbox, axis, userdef, undefined};
enum BoundarySide { left, right};
enum class SliceType {Cut, Average};

// Type of grid coarsening
enum GridCoarsening{disabled,
                    enabled,
                    dynamic}; ///< enabled = static coarsening (static is a reserved c++ keyword)

// Commonly used classes and functions
#include "global.hpp"
#include "error.hpp"
#include "macros.hpp"
#include "loop.hpp"
#include "reduce.hpp"
#include "arrays.hpp"

#endif  //  IDEFIX_HPP_
