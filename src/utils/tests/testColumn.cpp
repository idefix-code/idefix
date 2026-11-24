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

// only applies when being 3D + Spherical
#if DIMENSIONS == 3 && GEOMETRY == SPHERICAL
  // The test fill a 3D var with constant value and applies
  // Colum on it. Along theta it should have constant increase
  TEST(Column, ComputeColumn_4D_const_increase_along_theta) {
    //ini idefix app
    idfx::UnitTestApp app(8, 8, 8, 1.0);

    //extract one pointer for use
    DataBlock & data = *app.dataBlock;

    //spawn an array
    const int MYVAR = 0;
    const int VARS = 1;
    IdefixArray4D<real> vars = IdefixArray4D<real>("myvar", VARS, data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);

    //fill with const values
    ArrayUniformFill(vars, 3.0);

    //create column object & use
    Column column(JDIR, 1, &data);//teta
    column.ComputeColumn(vars, MYVAR);
    IdefixArray3D<real> result = column.GetColumn();

    //copy to ease checking
    IdefixHostArray3D<real> tmp("result_host", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);
    Kokkos::deep_copy(tmp, result);

    //should be const increase along theta except having 0 on ghosts
    for (int i = 2 ; i < data.np_tot[IDIR] - 2 ; i++) {
      for (int k = 2 ; k < data.np_tot[KDIR] - 2  ; k++) {
        double prev = 0.0;
        for (int j = 0 ; j < data.np_tot[JDIR] ; j++) {
          if (j < 2 || j >= 10) {
            ASSERT_EQ(0, tmp(k,j,i)) << "Ghost";
          } else if (j == 2) {
            prev = tmp(k,j,i);
            ASSERT_NE(0, tmp(k,j,i)) << "First non ghost";
          } else {
            ASSERT_DOUBLE_EQ(prev * (j-1), tmp(k,j,i)) << "Shound have constant increase along theta";
          }
        }
      }
    }
  }

  // The test fill a 3D var with constant value and applies
  // Colum on it. Along phi it should have constant increase
  TEST(Column, ComputeColumn_4D_const_increase_along_phi) {
    //ini idefix app
    idfx::UnitTestApp app(8, 8, 8, 1.0);

    //extract one pointer for use
    DataBlock & data = *app.dataBlock;

    //spawn an array
    const int MYVAR = 0;
    const int VARS = 1;
    IdefixArray4D<real> vars = IdefixArray4D<real>("myvar", VARS, data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);

    //fill with const values
    ArrayUniformFill(vars, 3.0);

    //create column object & use
    Column column(KDIR, 1, &data);//phi
    column.ComputeColumn(vars, MYVAR);
    IdefixArray3D<real> result = column.GetColumn();

    //copy on host to ease comare
    IdefixHostArray3D<real> tmp("result_host", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);
    Kokkos::deep_copy(tmp, result);

    //should be const increase along phi except having 0 on ghost
    for (int i = 2 ; i < data.np_tot[IDIR] - 2 ; i++) {
      for (int j = 2 ; j < data.np_tot[JDIR] - 2  ; j++) {
        double prev = 0.0;
        for (int k = 0 ; k < data.np_tot[KDIR] ; k++) {
          if (k < 2 || k >= 10) {
            ASSERT_EQ(0, tmp(k,j,i)) << "Ghost";
          } else if (k == 2) {
            prev = tmp(k,j,i);
            ASSERT_NE(0, tmp(k,j,i)) << "First non ghost";
          } else {
            ASSERT_DOUBLE_EQ(prev * (k-1), tmp(k,j,i)) << "Shound have constant increase along phi";
          }
        }
      }
    }
  }

  // The test fill a 3D var with constant value and applies
  // Colum on it. Along R it should have the pre-defined values
  // + being constant over theta & phi
  TEST(Column, ComputeColumn_4D_along_r_and_const_along_theta_phi) {
    //ini idefix app
    idfx::UnitTestApp app(8, 8, 8, 1.0);

    //extract one pointer for use
    DataBlock & data = *app.dataBlock;

    //spawn an array
    const int MYVAR = 0;
    const int VARS = 1;
    IdefixArray4D<real> vars = IdefixArray4D<real>("myvar", VARS, data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);

    //fill with const values
    ArrayUniformFill(vars, 3.0);

    //create column object & use
    Column column(IDIR, 1, &data);//r
    column.ComputeColumn(vars, MYVAR);
    IdefixArray3D<real> result = column.GetColumn();

    //copy on host to each checking
    IdefixHostArray3D<real> tmp("result_host", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);
    Kokkos::deep_copy(tmp, result);

    //should increase along R and be constant over phi/theta
    double expected_along_r[] = {
      0, 0,
      2.8, 5.7230769230769241, 8.6830769230769249, 11.658686679174485,
      14.642293236551534, 17.630528530669181, 20.621678973147059, 23.61478242142292,
      0, 0
    };
    for (int k = 2 ; k < data.np_tot[KDIR] - 2 ; k++)
      for (int j = 2 ; j < data.np_tot[JDIR] - 2  ; j++)
        for (int i = 0 ; i < data.np_tot[IDIR] ; i++)
          ASSERT_DOUBLE_EQ(expected_along_r[i], tmp(k,j,i)) << "Shound have constant along phi & theta";
  }

  // The test fill a 3D var with constant value and applies
  // Colum on it. Along R it should have the pre-defined values
  // + being constant over theta & phi
  TEST(Column, ComputeColumn_3D_along_r_and_const_along_theta_phi) {
    //ini idefix app
    idfx::UnitTestApp app(8, 8, 8, 1.0);

    //extract one pointer for use
    DataBlock & data = *app.dataBlock;

    //spawn an array
    IdefixArray3D<real> var = IdefixArray3D<real>("myvar", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);

    //fill with const values
    ArrayUniformFill(var, 3.0);

    //create column object & use
    Column column(IDIR, 1, &data);//r
    column.ComputeColumn(var);
    IdefixArray3D<real> result = column.GetColumn();

    //copy on host to each checking
    IdefixHostArray3D<real> tmp("result_host", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);
    Kokkos::deep_copy(tmp, result);

    //should increase along R and be constant over phi/theta
    double expected_along_r[] = {
      0, 0,
      2.8, 5.7230769230769241, 8.6830769230769249, 11.658686679174485,
      14.642293236551534, 17.630528530669181, 20.621678973147059, 23.61478242142292,
      0, 0
    };
    for (int k = 2 ; k < data.np_tot[KDIR] - 2 ; k++)
      for (int j = 2 ; j < data.np_tot[JDIR] - 2  ; j++)
        for (int i = 0 ; i < data.np_tot[IDIR] ; i++)
          ASSERT_DOUBLE_EQ(expected_along_r[i], tmp(k,j,i)) << "Shound have constant along phi & theta";
  }
#endif //DIMENSION == 3 && GEOMETRY == SPHERICAL
