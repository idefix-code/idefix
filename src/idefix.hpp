#ifndef IDEFIX_HPP
#define IDEFIX_HPP
#include <iostream>
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

using Device = Kokkos::DefaultExecutionSpace;
using Layout = Kokkos::LayoutRight;



//#define USE_DOUBLE

#include "real_types.h"
#include "loop.hpp"
#include "arrays.hpp"
#include "globals.hpp"
#include "input.hpp"
#include "grid.hpp"
//#include "configuration.hpp"
//#include "grid.hpp"


// Some macro definitions

#define     IDIR    0
#define     JDIR    1
#define     KDIR    2



#endif
