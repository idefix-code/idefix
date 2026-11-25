// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/tests/testArrays.hpp
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
#include "idefix.hpp"
#include "arrays.hpp"
#include "gtest-idfx.hpp"

TEST(Array, assert_1d) {
  //data array
  IdefixArray1D<double> a("a", 10);

  //ref array
  IdefixArray1D<double> b("b", 10);

  //make both 5
  idefix_for("Loop1", 0, 10, KOKKOS_LAMBDA (const int i) {
    a(i) = b(i) = 5;
  });

  //check all 5 static
  ASSERT_THAT_1D_ARRAY(a, testing::Eq, 5);

  //compare both array
  ASSERT_THAT_1D_ARRAYS(a, testing::Eq, b);

  //compare with near
  ASSERT_THAT_1D_ARRAYS_ARGS(a, testing::DoubleNear, b, 0.001);
}

TEST(Array, assert_1d_is_nan) {
  //data array
  IdefixArray1D<double> a("a", 10);

  //ref array
  IdefixArray1D<double> b("b", 10);

  //make both 5
  idefix_for("Loop1", 0, 10, KOKKOS_LAMBDA (const int i) {
    a(i) = b(i) = NAN;
  });

  //check all 5 static
  ASSERT_THAT_1D_ARRAY(a, testing::IsNan);
}

TEST(Array, assert_2d) {
  //data array
  IdefixArray2D<double> a("a", 10, 20);

  //ref array
  IdefixArray2D<double> b("b", 10, 20);

  //make both 5
  idefix_for("Loop1", 0, 10, 0, 20, KOKKOS_LAMBDA (const int j, const int i) {
    a(j, i) = b(j, i) = 5;
  });

  //check all 5 static
  ASSERT_THAT_2D_ARRAY(a, testing::Eq, 5);

  //compare both array
  ASSERT_THAT_2D_ARRAYS(a, testing::Eq, b);

  //compare with near
  ASSERT_THAT_2D_ARRAYS_ARGS(a, testing::DoubleNear, b, 0.001);
}

TEST(Array, assert_3d) {
  //data array
  IdefixArray3D<double> a("a", 10, 20, 30);

  //ref array
  IdefixArray3D<double> b("b", 10, 20, 30);

  //make both 5
  idefix_for("Loop1", 0, 10, 0, 20, 0, 30, KOKKOS_LAMBDA (const int k, const int j, const int i) {
    a(k, j, i) = b(k, j, i) = 5;
  });

  //check all 5 static
  ASSERT_THAT_3D_ARRAY(a, testing::Eq, 5);

  //compare both array
  ASSERT_THAT_3D_ARRAYS(a, testing::Eq, b);

  //compare with near
  ASSERT_THAT_3D_ARRAYS_ARGS(a, testing::DoubleNear, b, 0.001);
}

TEST(Array, assert_4d) {
  //data array
  IdefixArray4D<double> a("a", 10, 20, 30, 40);

  //ref array
  IdefixArray4D<double> b("b", 10, 20, 30, 40);

  //make both 5
  idefix_for("Loop1", 0, 10, 0, 20, 0, 30, 0, 40,
    KOKKOS_LAMBDA (const int l, const int k, const int j, const int i) {
    a(l, k, j, i) = b(l, k, j, i) = 5;
  });

  //check all 5 static
  ASSERT_THAT_4D_ARRAY(a, testing::Eq, 5);

  //compare both array
  ASSERT_THAT_4D_ARRAYS(a, testing::Eq, b);

  //compare with near
  ASSERT_THAT_4D_ARRAYS_ARGS(a, testing::DoubleNear, b, 0.001);
}
