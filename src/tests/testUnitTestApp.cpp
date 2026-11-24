// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/tests/testUnitTestApp.hpp
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
#include "gtest-idfx.hpp"
#include "unitTestApp.hpp"

TEST(UnitTestApp, constructor) {
  //ini idefix app
  idfx::UnitTestApp app(8,8,8);

  //extract one pointer for use
  DataBlock & data = *app.dataBlock;

  //checks about mesh with ghost cells
  ASSERT_EQ(data.np_tot[IDIR], 12);
  ASSERT_EQ(data.np_tot[JDIR], 12);
  ASSERT_EQ(data.np_tot[KDIR], 12);

  //checks about internal mesh
  ASSERT_EQ(data.np_int[IDIR], 8);
  ASSERT_EQ(data.np_int[JDIR], 8);
  ASSERT_EQ(data.np_int[KDIR], 8);

  //checks about cell size
  ASSERT_EQ(data.dx[IDIR][0], 1.0);
  ASSERT_EQ(data.dx[JDIR][0], 1.0);
  ASSERT_EQ(data.dx[KDIR][0], 1.0);
}
