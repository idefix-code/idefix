#ifndef IDEFIX_HPP
#define IDEFIX_HPP
#include <iostream>
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

using Device = Kokkos::DefaultExecutionSpace;
using Layout = Kokkos::LayoutRight;

// Define the kind of loop we want (see loop.hpp for details)
//#define  INDEX_LOOP
//#define  MDRANGE_LOOP
//#define  SIMD_LOOP
//#define  TP_INNERX_LOOP
//#define  TPTTRTVR_LOOP

// Hopefully a master switch which detects which loop is needed on the architecture

#ifdef KOKKOS_ENABLE_CUDA
#define INDEX_LOOP
#else
#define TP_INNERX_LOOP
#endif

#define USE_DOUBLE

#define     YES     255
#define     NO      0

// Basic configuration
#include "real_types.h"
#include "definitions.hpp"


// Shortcuts for fields used in the code

#define  RHO 0
#define  MX1 1
#define  MX2 (COMPONENTS >= 2 ? 2: 255)
#define  MX3 (COMPONENTS == 3 ? 3: 254)
#if MHD == YES
#define  BX1 (COMPONENTS + 1)
#define  BX2 (COMPONENTS >= 2 ? (BX1+1): 252)
#define  BX3 (COMPONENTS == 3 ? (BX1+2): 251)
#else
#define  BX1 253
#define  BX2 252
#define  BX3 251
#endif

#if HAVE_ENERGY
#if MHD == YES
  #define ENG  (2*COMPONENTS + 1)
#else
  #define ENG  (COMPONENTS + 1)
#endif
  #define PRS  ENG
#endif

#define VX1   MX1
#define VX2   MX2
#define VX3   MX3

#if MHD == YES
  #define NFLX (1 + 2*COMPONENTS + HAVE_ENERGY)
#else
  #define NFLX (1 + COMPONENTS + HAVE_ENERGY)
#endif

// Face-centered variables
#define BX1s  0
#define BX2s  1
#define BX3s  2


// Some macro definitions

/* ---- Geometry Labels ( > 0) ----  */

#define CARTESIAN    1
#define CYLINDRICAL  2
#define POLAR        3
#define SPHERICAL    4

#define     IDIR    0
#define     JDIR    1
#define     KDIR    2



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

/*---- EMFs -----*/
#define ARITHMETIC   1
#define UCT0         2
#define UCT_CONTACT  3

#ifndef EMF_AVERAGE
    #define EMF_AVERAGE     ARITHMETIC
#endif


// Types of boundary which can be treated
enum BoundaryType { internal, periodic, outflow, shearingbox, userdef};
enum BoundarySide { left, right};



#include "error.hpp"
#include "macros.hpp"
#include "loop.hpp"
#include "arrays.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "gridHost.hpp"
#include "dataBlock.hpp"
#include "dataBlockHost.hpp"
#include "test.hpp"
#include "hydro.hpp"
#include "timeIntegrator.hpp"
#include "setup.hpp"
#include "outputVtk.hpp"

#ifndef MHD
#error MHD flag should be set to yes or no in definitions.hpp
#endif


#endif
