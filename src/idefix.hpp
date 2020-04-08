#ifndef IDEFIX_HPP
#define IDEFIX_HPP
#include <iostream>
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

using Device = Kokkos::DefaultExecutionSpace;
using Layout = Kokkos::LayoutRight;

// Define the kind of loop we want (see loop.hpp for details)
//#define  INDEX_LOOP
#define  MDRANGE_LOOP
//#define  SIMD_LOOP
//#define  TP_INNERX_LOOP
//#define  TPTTRTVR_LOOP

#define USE_DOUBLE






//#include "configuration.hpp"
//#include "grid.hpp"


// Some macro definitions

/* ---- Geometry Labels ( > 0) ----  */

#define CARTESIAN    1
#define CYLINDRICAL  2
#define POLAR        3
#define SPHERICAL    4

#define     IDIR    0
#define     JDIR    1
#define     KDIR    2

#define     YES     255
#define     NO      0

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


// Types of boundary which can be treated
enum BoundaryType { internal, periodic, outflow, userdef};
enum BoundarySide { left, right};


#include "real_types.h"
#include "definitions.hpp"
#include "error.hpp"
#include "macros.hpp"
#include "mod_defs.hpp"
#include "loop.hpp"
#include "arrays.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "gridHost.hpp"
#if MHD == YES
#include "electromotiveforce.hpp"
#endif
#include "dataBlock.hpp"
#include "dataBlockHost.hpp"
#include "test.hpp"
#include "setup.hpp"
#include "physics.hpp"
#include "timeIntegrator.hpp"
#include "outputVtk.hpp"



#endif
