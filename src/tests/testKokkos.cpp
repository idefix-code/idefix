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
  Kokkos::View<double*, Kokkos::DefaultExecutionSpace> a("a", 10);

  //ref array
  Kokkos::View<double*, Kokkos::DefaultExecutionSpace> b("a", 10);

  //make both 5
  for (size_t i = 0 ; i < 10 ; i++) {
    a(i) = b(i) = 5;
  }

  //check all 5 static
  ASSERT_THAT_1D_ARRAY(a, testing::Eq, 5);

  //compare both array
  ASSERT_THAT_1D_ARRAYS(a, testing::Eq, b);

  //compare with near
  ASSERT_THAT_1D_ARRAYS_ARGS(a, testing::DoubleNear, b, 0.001);
}

TEST(Kokkos, assert_view_1d_is_nan) {
  //data array
  Kokkos::View<double*, Kokkos::DefaultExecutionSpace> a("a", 10);

  //ref array
  Kokkos::View<double*, Kokkos::DefaultExecutionSpace> b("a", 10);

  //make both 5
  for (size_t i = 0 ; i < 10 ; i++) {
    a(i) = b(i) = NAN;
  }

  //check all 5 static
  ASSERT_THAT_1D_ARRAY(a, testing::IsNan);
}

TEST(Kokkos, assert_view_2d) {
  //data array
  Kokkos::View<double**, Kokkos::DefaultExecutionSpace> a("a", 10, 10);

  //ref array
  Kokkos::View<double**, Kokkos::DefaultExecutionSpace> b("a", 10, 10);

  //make both 5
  for (size_t j = 0 ; j < 10 ; j++) {
    for (size_t i = 0 ; i < 10 ; i++) {
      a(i, j) = b(i, j) = 5;
    }
  }

  //check all 5 static
  ASSERT_THAT_2D_ARRAY(a, testing::Eq, 5);

  //compare both array
  ASSERT_THAT_2D_ARRAYS(a, testing::Eq, b);

  //compare with near
  ASSERT_THAT_2D_ARRAYS_ARGS(a, testing::DoubleNear, b, 0.001);
}

TEST(Kokkos, assert_view_3d) {
  //data array
  Kokkos::View<double***, Kokkos::DefaultExecutionSpace> a("a", 10, 10, 10);

  //ref array
  Kokkos::View<double***, Kokkos::DefaultExecutionSpace> b("a", 10, 10, 10);

  //make both 5
  for (size_t k = 0 ; k < 10 ; k++) {
    for (size_t j = 0 ; j < 10 ; j++) {
      for (size_t i = 0 ; i < 10 ; i++) {
        a(i, j, k) = b(i, j, k) = 5;
      }
    }
  }

  //check all 5 static
  ASSERT_THAT_3D_ARRAY(a, testing::Eq, 5);

  //compare both array
  ASSERT_THAT_3D_ARRAYS(a, testing::Eq, b);

  //compare with near
  ASSERT_THAT_3D_ARRAYS_ARGS(a, testing::DoubleNear, b, 0.001);
}

TEST(Array, assert_view_4d) {
  //data array
  Kokkos::View<double****, Kokkos::DefaultExecutionSpace> a("a", 10, 10, 10, 10);

  //ref array
  Kokkos::View<double****, Kokkos::DefaultExecutionSpace> b("b", 10, 10, 10, 10);

  //make both 5
  for (size_t l = 0 ; l < 10 ; l++) {
    for (size_t k = 0 ; k < 10 ; k++) {
      for (size_t j = 0 ; j < 10 ; j++) {
        for (size_t i = 0 ; i < 10 ; i++) {
          a(i, j, k, l) = b(i, j, k, l) = 5;
        }
      }
    }
  }

  //check all 5 static
  ASSERT_THAT_4D_ARRAY(a, testing::Eq, 5);

  //compare both array
  ASSERT_THAT_4D_ARRAYS(a, testing::Eq, b);

  //compare with near
  ASSERT_THAT_4D_ARRAYS_ARGS(a, testing::DoubleNear, b, 0.001);
}
