#ifndef IDEFIX_HPP
#define IDEFIX_HPP
#include <iostream>
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

using Device = Kokkos::DefaultExecutionSpace;
using Layout = Kokkos::LayoutRight;

// Define the kind of loop we want (see loop.hpp for details)
#define  INDEX_LOOP

// #define USE_DOUBLE

#include "real_types.h"
#include "definitions.hpp"
#include "mod_defs.hpp"
#include "loop.hpp"
#include "arrays.hpp"
#include "globals.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "gridHost.hpp"
#include "data.hpp"
#include "dataHost.hpp"
#include "test.hpp"




//#include "configuration.hpp"
//#include "grid.hpp"


// Some macro definitions


#define     IDIR    0
#define     JDIR    1
#define     KDIR    2


#define NVAR    (NFLX)




#endif
