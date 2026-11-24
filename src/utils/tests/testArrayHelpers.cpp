// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/utils/tests/testArrayHelpers.hpp
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
#include "idefix.hpp"
#include "arrays.hpp"
#include "arrayHelpers.hpp"
#include "gtest-idfx.hpp"

TEST(ArrayHelpers, ArrayUniformFill_1D) {
  //array
  IdefixArray1D<real> array("array", 10);

  //fill with uniform
  ArrayUniformFill(array, 5.0);

  //check
  ASSERT_THAT_1D_ARRAY(array, testing::Eq, 5.0);
}

TEST(ArrayHelpers, ArrayUniformFill_2D) {
  //array
  IdefixArray2D<real> array("array", 10, 10);

  //fill with uniform
  ArrayUniformFill(array, 5.0);

  //check
  ASSERT_THAT_2D_ARRAY(array, testing::Eq, 5.0);
}

TEST(ArrayHelpers, ArrayUniformFill_3D) {
  //array
  IdefixArray3D<real> array("array", 10, 10, 10);

  //fill with uniform
  ArrayUniformFill(array, 5.0);

  //check
  ASSERT_THAT_3D_ARRAY(array, testing::Eq, 5.0);
}

TEST(ArrayHelpers, ArrayUniformFill_4D) {
  //array
  IdefixArray4D<real> array("array", 10, 10, 10, 10);

  //fill with uniform
  ArrayUniformFill(array, 5.0);

  //check
  ASSERT_THAT_4D_ARRAY(array, testing::Eq, 5.0);
}
