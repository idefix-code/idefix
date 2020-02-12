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
#define  SIMD_LOOP
//#define  TP_INNERX_LOOP
//#define  TPTTRTVR_LOOP

#define USE_DOUBLE

#include "real_types.h"
#include "definitions.hpp"
#include "mod_defs.hpp"
#include "loop.hpp"
#include "arrays.hpp"
#include "globals.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "gridHost.hpp"
#include "dataBlock.hpp"
#include "dataBlockHost.hpp"
#include "test.hpp"
#include "timeIntegrator.hpp"




//#include "configuration.hpp"
//#include "grid.hpp"


// Some macro definitions


#define     IDIR    0
#define     JDIR    1
#define     KDIR    2


#define NVAR    (NFLX)




#endif
