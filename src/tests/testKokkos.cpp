// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/tests/testKokkos.hpp
//
// Last modified : 10/2026
//
// Copyright(C) by :
// - Sébastien Valat <sebastien.valat@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// headers
#include <vector>
#include <string>
#include <Kokkos_Core.hpp>
#include "gtest-idfx.hpp"

TEST(Kokkos, assert_view_1d) {
  //data array
  Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> a("a", 10);

  //ref array
  Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> b("a", 10);

  //make both 5
  Kokkos::parallel_for("Loop1", 10, KOKKOS_LAMBDA (const int i) {
    a(i) = b(i) = 5;
  });

  //check all 5 static
  ASSERT_THAT_1D_ARRAY(a, testing::Eq, 5);

  //compare both array
  ASSERT_THAT_1D_ARRAYS(a, testing::Eq, b);

  //compare with near
  ASSERT_THAT_1D_ARRAYS_ARGS(a, testing::DoubleNear, b, 0.001);
}

TEST(Kokkos, assert_view_1d_is_nan) {
  //data array
  Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> a("a", 10);

  //ref array
  Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> b("a", 10);

  //make both 5
  Kokkos::parallel_for("Loop1", 10, KOKKOS_LAMBDA (const int i) {
    a(i) = b(i) = NAN;
  });

  //check all 5 static
  ASSERT_THAT_1D_ARRAY(a, testing::IsNan);
}
