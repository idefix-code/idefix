// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/utils/tests/testColumn.hpp
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
#include "column.hpp"
#include "arrayHelpers.hpp"

TEST(UtIdefixEnv, ComputeColumn_X1_Right) {
  //ini idefix app
  idfx::UnitTestApp app(8, 8, 8, 1.0);

  //extract one pointer for use
  DataBlock & data = *app.dataBlock;

  //spawn an array
  const int MYVAR = 0;
  const int VARS = 1;
  IdefixArray4D<real> vars = IdefixArray4D<real>("myvar", VARS, data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);

  //fill
  ArrayUniformFill(vars, 1.0);

  //build expected array
  IdefixArray3D<real> expected = IdefixArray3D<real>("expected", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);

  //create column object & use
  Column column(IDIR, -1, &data);
  column.ComputeColumn(vars, MYVAR);
  IdefixArray3D<real> result = column.GetColumn();

  //check
  ASSERT_THAT_3D_ARRAYS(result, testing::Eq, expected);
}
