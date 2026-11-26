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
    Kokkos::deep_copy(vars, 3.0);

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
            const real r = data.x[0](i);
            const real dtheta = 3.1415926/8.0;
            //const real dV = r * r * sin(theta) * dr * dtheta * dphi;
            //const real avgA = r * sin(theta) * dr * dphi;
            //const real expected = 3.0 * dV / avgA;
            const real expected = 3.0 * r * dtheta;
            ASSERT_NEAR(expected, tmp(k,j,i), 0.15) << "First non ghost with single contribution";

            //remind the real value to check the constant increase
            prev = tmp(k,j,i);
          } else {
            ASSERT_NEAR(prev * (j-1), tmp(k,j,i),10e-14) << "Shound have constant increase along theta";
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
    DataBlockHost & data = *app.dataBlockHost;

    //spawn an array
    const int MYVAR = 0;
    const int VARS = 1;
    IdefixArray4D<real> vars = IdefixArray4D<real>("myvar", VARS, data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);

    //fill with const values
    Kokkos::deep_copy(vars, 3.0);

    //create column object & use
    Column column(KDIR, 1, app.dataBlock.get());//phi
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
            const real r = data.x[0](i);
            const real theta = data.x[1](j);
            const real dphi = 2*3.1415926/8.0;
            //const real dV = r * r * sin(theta) * dr * dtheta * dphi;
            //const real avgA = r * dr * dtheta;
            //const real expected = 3.0 * dV / avgA;
            const real expected = 3.0 * r * sin(theta) * dphi;
            ASSERT_NEAR(expected, tmp(k,j,i), 0.11) << "First non ghost with single contribution";
            //remind the real value to check the constant increase
            prev = tmp(k,j,i);
          } else {
            ASSERT_NEAR(prev * (k-1), tmp(k,j,i), 10e-14) << "Shound have constant increase along phi";
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
    DataBlockHost & data = *app.dataBlockHost;

    //spawn an array
    const int MYVAR = 0;
    const int VARS = 1;
    IdefixArray4D<real> vars = IdefixArray4D<real>("myvar", VARS, data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);

    //fill with const values
    Kokkos::deep_copy(vars, 3.0);

    //create column object & use
    Column column(IDIR, 1, app.dataBlock.get());//r
    column.ComputeColumn(vars, MYVAR);
    IdefixArray3D<real> result = column.GetColumn();

    //copy on host to each checking
    IdefixHostArray3D<real> tmp("result_host", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);
    Kokkos::deep_copy(tmp, result);

    //should increase along R and be constant over phi/theta
    for (int i = 0 ; i < data.np_tot[IDIR] ; i++) {
      real prev = 0.0;
      for (int j = 2 ; j < data.np_tot[JDIR] - 2  ; j++) {
        for (int k = 2 ; k < data.np_tot[KDIR] - 2 ; k++) {
          if (i < 2 || i >= 10) {
            ASSERT_EQ(0, tmp(k,j,i)) << "Ghost";
          } else if (prev == 0.0) {
            const real dr = 1.0;
            //const real dV = r * r * sin(theta) * dr * dtheta * dphi;
            //const real avgA = r * r * sin(theta) * dtheta * dphi;
            //const real expected = (i-1) * dV / avgA
            const real expected = 3.0 * dr * (i-1);
            ASSERT_NEAR(expected, tmp(k,j,i), 0.4) << "First non ghost with single contribution";
            //remind the real value to check the constant value over theta phi
            prev = tmp(k,j,i);
          } else {
            ASSERT_NEAR(prev, tmp(k,j,i), 10e-14) << "Shound have constant increase along theta phi";
          }
        }
      }
    }
  }

  // The test fill a 3D var with constant value and applies
  // Colum on it. Along R it should have the pre-defined values
  // + being constant over theta & phi
  TEST(Column, ComputeColumn_3D_along_r_and_const_along_theta_phi) {
    //ini idefix app
    idfx::UnitTestApp app(8, 8, 8, 1.0);

    //extract one pointer for use
    DataBlockHost & data = *app.dataBlockHost;

    //spawn an array
    IdefixArray3D<real> var = IdefixArray3D<real>("myvar", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);

    //fill with const values
    Kokkos::deep_copy(var, 3.0);

    //create column object & use
    Column column(IDIR, 1, app.dataBlock.get());//r
    column.ComputeColumn(var);
    IdefixArray3D<real> result = column.GetColumn();

    //copy on host to each checking
    IdefixHostArray3D<real> tmp("result_host", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);
    Kokkos::deep_copy(tmp, result);

    //should increase along R and be constant over phi/theta
    for (int i = 0 ; i < data.np_tot[IDIR] ; i++) {
      real prev = 0.0;
      for (int j = 2 ; j < data.np_tot[JDIR] - 2  ; j++) {
        for (int k = 2 ; k < data.np_tot[KDIR] - 2 ; k++) {
          if (i < 2 || i >= 10) {
            ASSERT_EQ(0, tmp(k,j,i)) << "Ghost";
          } else if (prev == 0.0) {
            const real dr = 1.0;
            //const real dV = r * r * sin(theta) * dr * dtheta * dphi;
            //const real avgA = r * r * sin(theta) * dtheta * dphi;
            //const real expected = (i-1) * dV / avgA
            const real expected = 3.0 * dr * (i-1);
            ASSERT_NEAR(expected, tmp(k,j,i), 0.4) << "First non ghost with single contribution";
            //remind the real value to check the constant value over theta phi
            prev = tmp(k,j,i);
          } else {
            ASSERT_NEAR(prev, tmp(k,j,i), 10e-14) << "Should have constant value over phi / theta";
          }
        }
      }
    }
  }
#endif //DIMENSION == 3 && GEOMETRY == SPHERICAL
