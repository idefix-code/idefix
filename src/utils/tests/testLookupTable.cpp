// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/utils/testLookupTable.hpp
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
#include "gtest-idfx.hpp"
#include "lookupTable.hpp"

// Custom transformation (and its inverse) to check that the functions used for the interpolation
// in function space can be changed when the lookup table is created
struct MyLog10 {
  KOKKOS_INLINE_FUNCTION real operator() (const real x) const {
    return(log10(x));
  }
};

struct MyPow10 {
  KOKKOS_INLINE_FUNCTION real operator() (const real x) const {
    return(pow(10.0,x));
  }
};

// Expected bracketing neighbours of the 2D CSV table (toto.csv) around (x=2.1, y=3.5): x brackets
// [2,3], y brackets [3,4], and the data on the four surrounding vertices is 5,6,6,7. This is
// checked from several angles across the file below (Get vs GetNeighbours vs GetNeighboursIndx,
// host vs device, cached vs uncached search), so it is defined once here instead of being
// retyped as magic numbers in every block.
const real kXN2D[4]      = {2.0, 3.0, 3.0, 4.0};
const real kDataN2D[4]   = {5.0, 6.0, 6.0, 7.0};
const int  kIdx2D[2]     = {0, 1};
const int  kDataIdx2D[4] = {1, 4, 2, 5};
const real kValue2D      = 5.6;   // csv.Get({2.1, 3.5})

//Testing 2D CSV file on device
TEST(lookupTable, 2d_csv_file_on_device) {
  IdefixArray1D<real> arr = IdefixArray1D<real>("Test",1);
  IdefixArray1D<real>::host_mirror_type arrHost = Kokkos::create_mirror_view(arr);

  LookupTable<2> csv(DATA_FILE("toto.csv"),',');

  idefix_for("loop",0, 1, KOKKOS_LAMBDA (int i) {
    real x[2];
    x[0] = 2.1;
    x[1] = 3.5;
    arr(i) = csv.Get(x);
  });

  Kokkos::deep_copy(arrHost , arr);

  ASSERT_NEAR(arrHost(0), kValue2D, 1e-13) <<  "2D CSV, device";
}

//Testing 2D CSV file on Host
TEST(lookupTable, 2d_csv_file_on_host) {
    LookupTable<2> csv(DATA_FILE("toto.csv"),',');
    real x[2];
    x[0] = 2.1;
    x[1] = 3.5;
    real result = csv.GetHost(x);
    ASSERT_NEAR(result, kValue2D, 1e-13) << "2D CSV, host";
}

//Testing 1D CSV file on device.
TEST(lookupTable, 1d_csv_file_on_device) {
  IdefixArray1D<real> arr = IdefixArray1D<real>("Test",1);
  IdefixArray1D<real>::host_mirror_type arrHost = Kokkos::create_mirror_view(arr);

  // Read 1D CSV File
  LookupTable<1> csv1D(DATA_FILE("toto1D.csv"),',');

  idefix_for("loop",0, 1, KOKKOS_LAMBDA (int i) {
    real x[2];
    x[0] = 2.1;
    x[1] = 3.5;
    arr(i) = csv1D.Get(x);
  });

  Kokkos::deep_copy(arrHost , arr);
  ASSERT_EQ(arrHost(0), real(4.2)) << "1D CSV, device";
}

//Testing 1D CSV file on host.
TEST(lookupTable, 1d_csv_file_on_host) {
  LookupTable<1> csv1D(DATA_FILE("toto1D.csv"),',');
  real x[2];
  x[0] = 2.1;
  x[1] = 3.5;
  real result = csv1D.GetHost(x);
  ASSERT_EQ(result, real(4.2)) << "1D CSV, host";
}

//Testing 1D CSV file read as columns on device.
TEST(LookupTable, 1d_csv_read_as_col_on_device) {
  // Read the same 1D table, but stored as columns of the CSV file
  LookupTable<1> csv1Dcolumn(DATA_FILE("toto1Dcolumn.csv"),',', true, true);
  IdefixArray1D<real> arr = IdefixArray1D<real>("Test",1);
  IdefixArray1D<real>::host_mirror_type arrHost = Kokkos::create_mirror_view(arr);

  idefix_for("loop",0, 1, KOKKOS_LAMBDA (int i) {
    real x[1];
    x[0] = 2.1;
    arr(i) = csv1Dcolumn.Get(x);
  });

  Kokkos::deep_copy(arrHost , arr);
  ASSERT_EQ(arrHost(0), real(4.2)) << "1D CSV read as columns, device";
}

//Testing 1D CSV file read as columns on Host.
TEST(LookupTable, 1d_csv_read_as_col_on_host) {
  // Read the same 1D table, but stored as columns of the CSV file
  LookupTable<1> csv1Dcolumn(DATA_FILE("toto1Dcolumn.csv"),',', true, true);

  real x[1];
  x[0] = 2.1;
  real result = csv1Dcolumn.GetHost(x);
  ASSERT_EQ(result, real(4.2)) << "1D CSV read as columns, host";
}

//Testing interpolation in function space on device.
TEST(LookupTable, interpol_func_space_on_device) {
  IdefixArray1D<real> arr = IdefixArray1D<real>("Test",1);
  IdefixArray1D<real>::host_mirror_type arrHost = Kokkos::create_mirror_view(arr);

  // data = x^2 sampled in x=1,2,4. A linear interpolation in x=3 gives 10, while an
  // interpolation performed on log(x) and log(data) gives the exact result 9
  LookupTable<1> csv1Dlin(DATA_FILE("toto1Dlog.csv"),',', true, true);
  LookupTable<1> csv1Dlog(DATA_FILE("toto1Dlog.csv"),',', true, true, true);
  // the transformation and its inverse can also be chosen when the table is created
  LookupTable<1, MyLog10, MyPow10> csv1Dlog10(DATA_FILE("toto1Dlog.csv"),',', true, true, true);

  idefix_for("loop",0, 1, KOKKOS_LAMBDA (int i) {
    real x[1];
    x[0] = 3.0;
    arr(i) = csv1Dlog.Get(x);
  });
  Kokkos::deep_copy(arrHost , arr);
  ASSERT_NEAR(arrHost(0), 9.0, 1e-13) << "function-space interpolation (default log), device";

  idefix_for("loop",0, 1, KOKKOS_LAMBDA (int i) {
    real x[1];
    x[0] = 3.0;
    arr(i) = csv1Dlog10.Get(x);
  });
  Kokkos::deep_copy(arrHost , arr);
  ASSERT_NEAR(arrHost(0), 9.0, 1e-13) << "function-space interpolation (custom log10), device";
}

//Testing interpolation in function space on Host.
TEST(LookupTable, interpol_func_space_on_host) {
  // data = x^2 sampled in x=1,2,4. A linear interpolation in x=3 gives 10, while an
  // interpolation performed on log(x) and log(data) gives the exact result 9
  LookupTable<1> csv1Dlin(DATA_FILE("toto1Dlog.csv"),',', true, true);
  LookupTable<1> csv1Dlog(DATA_FILE("toto1Dlog.csv"),',', true, true, true);
  // the transformation and its inverse can also be chosen when the table is created
  LookupTable<1, MyLog10, MyPow10> csv1Dlog10(DATA_FILE("toto1Dlog.csv"),',', true, true, true);

  real x[1];
  x[0] = 3.0;
  real result = csv1Dlog.GetHost(x);
  ASSERT_NEAR(result, 9.0, 1e-13) << "function-space interpolation (default log), host";

  result = csv1Dlog10.GetHost(x);
  ASSERT_NEAR(result, 9.0, 1e-13) << "function-space interpolation (custom log10), host";

  // the same table, interpolated linearly (default behaviour)
  result = csv1Dlin.GetHost(x);
  ASSERT_NEAR(result, 10.0, 1e-13) << "plain linear interpolation, host";
}

//Testing the neighbours used for the interpolation on Host.
TEST(LookupTable, neighbours_used_for_interpol_on_host) {
  // data = x^2 sampled in x=1,2,4. A linear interpolation in x=3 gives 10, while an
  // interpolation performed on log(x) and log(data) gives the exact result 9
  LookupTable<1> csv1Dlog(DATA_FILE("toto1Dlog.csv"),',', true, true, true);
  LookupTable<1> csv1D(DATA_FILE("toto1D.csv"),',');
  LookupTable<2> csv(DATA_FILE("toto.csv"),',');

  // 1D table (x=1,2,3 and data=2,4,6) interpolated in x=2.1
  real xq[1];
  xq[0] = 2.1;
  real xN[2];
  real dataN[2];
  int idx[1];
  int dataIdx[2];
  csv1D.GetNeighboursHost(xq, xN, dataN);
  csv1D.GetNeighboursIndxHost(xq, idx, dataIdx);
  //idfx::cout << "1D: x=" << xN[0] << "," << xN[1] << " data=" << dataN[0] << "," << dataN[1]
  //            << " idx=" << idx[0] << " dataIdx=" << dataIdx[0] << "," << dataIdx[1]
  //            << std::endl;
  const real expectedXN[2] = {2.0, 3.0};
  const real expectedDataN[2] = {4.0, 6.0};
  const int expectedDataIdx[2] = {1, 2};
  ASSERT_C_ARRAY_EQ(xN, expectedXN, 2, "1D neighbours, coordinates");
  ASSERT_C_ARRAY_EQ(dataN, expectedDataN, 2, "1D neighbours, data");
  ASSERT_EQ(idx[0], 1) << "1D neighbours, idx";
  ASSERT_C_ARRAY_EQ(dataIdx, expectedDataIdx, 2, "1D neighbours, dataIdx");

  // 2D table interpolated in (2.1,3.5). The vertices are ordered so that the bit n of the
  // vertex index tells whether we are on the right (1) or on the left (0) of the dimension n
  real xq2[2];
  xq2[0] = 2.1;
  xq2[1] = 3.5;
  real xN2[4];
  real dataN2[4];
  int idx2[2];
  int dataIdx2[4];
  csv.GetNeighboursHost(xq2, xN2, dataN2);
  csv.GetNeighboursIndxHost(xq2, idx2, dataIdx2);
  //idfx::cout << "2D: x=" << xN2[0] << "," << xN2[1] << " y=" << xN2[2] << "," << xN2[3]
  //            << " data=" << dataN2[0] << "," << dataN2[1] << "," << dataN2[2] << ","
  //            << dataN2[3] << " idx=" << idx2[0] << "," << idx2[1] << std::endl;
  ASSERT_C_ARRAY_EQ(xN2, kXN2D, 4, "2D neighbours, coordinates");
  ASSERT_C_ARRAY_EQ(dataN2, kDataN2D, 4, "2D neighbours, data");
  ASSERT_C_ARRAY_EQ(idx2, kIdx2D, 2, "2D neighbours, idx");
  ASSERT_C_ARRAY_EQ(dataIdx2, kDataIdx2D, 4, "2D neighbours, dataIdx");

  // The neighbours should reproduce the value returned by GetHost
  real dx = (xq2[0]-xN2[0])/(xN2[1]-xN2[0]);
  real dy = (xq2[1]-xN2[2])/(xN2[3]-xN2[2]);
  real interpolated = (1-dx)*(1-dy)*dataN2[0] + dx*(1-dy)*dataN2[1]
                    + (1-dx)*dy*dataN2[2] + dx*dy*dataN2[3];
  //idfx::cout << "2D: interpolation from the neighbours=" << interpolated << std::endl;
  ASSERT_NEAR(interpolated, csv.GetHost(xq2), 1e-13)
    << "2D neighbours, manual interpolation vs GetHost";

  // When the table is interpolated in function space, the neighbours are returned in the
  // original space of the table (x=1,2,4 and data=1,4,16 here)
  real xq3[1];
  xq3[0] = 3.0;
  real xN3[2];
  real dataN3[2];
  int idx3[1];
  int dataIdx3[2];
  csv1Dlog.GetNeighboursHost(xq3, xN3, dataN3);
  csv1Dlog.GetNeighboursIndxHost(xq3, idx3, dataIdx3);
  //idfx::cout << "1D (function space): x=" << xN3[0] << "," << xN3[1]
  //            << " data=" << dataN3[0] << "," << dataN3[1] << " idx=" << idx3[0] << std::endl;
  ASSERT_NEAR(xN3[0], 2.0, 1e-13) << "1D function-space neighbours, x lower bound";
  ASSERT_NEAR(xN3[1], 4.0, 1e-13) << "1D function-space neighbours, x upper bound";
  ASSERT_NEAR(dataN3[0], 4.0, 1e-13) << "1D function-space neighbours, data lower bound";
  ASSERT_NEAR(dataN3[1], 16.0, 1e-13) << "1D function-space neighbours, data upper bound";
  ASSERT_EQ(idx3[0], 1) << "1D function-space neighbours, idx";
  const int expectedDataIdx3[2] = {1, 2};
  ASSERT_C_ARRAY_EQ(dataIdx3, expectedDataIdx3, 2, "1D function-space neighbours, dataIdx");
}

//Testing the neighbours used for the interpolation on device.");
TEST(LookupTable, neighbours_used_for_interpol_on_device) {
  LookupTable<2> csv(DATA_FILE("toto.csv"),',');

  // data = x^2 sampled in x=1,2,4. A linear interpolation in x=3 gives 10, while an
  // interpolation performed on log(x) and log(data) gives the exact result 9
  LookupTable<1> csv1Dlin(DATA_FILE("toto1Dlog.csv"),',', true, true);
  LookupTable<1> csv1Dlog(DATA_FILE("toto1Dlog.csv"),',', true, true, true);

  IdefixArray1D<real> xNdev = IdefixArray1D<real>("xN",4);
  IdefixArray1D<real> dataNdev = IdefixArray1D<real>("dataN",4);
  IdefixArray1D<int> idxDev = IdefixArray1D<int>("idx",2);
  IdefixArray1D<int> dataIdxDev = IdefixArray1D<int>("dataIdx",4);

  idefix_for("neighbours",0, 1, KOKKOS_LAMBDA (int i) {
    real xq[2];
    xq[0] = 2.1;
    xq[1] = 3.5;
    real xN[4];
    real dataN[4];
    int idx[2];
    int dataIdx[4];
    csv.GetNeighbours(xq, xN, dataN);
    csv.GetNeighboursIndx(xq, idx, dataIdx);
    for(int n = 0 ; n < 4 ; n++) {
      xNdev(n) = xN[n];
      dataNdev(n) = dataN[n];
      dataIdxDev(n) = dataIdx[n];
    }
    idxDev(0) = idx[0];
    idxDev(1) = idx[1];
  });

  auto xNHost = Kokkos::create_mirror_view(xNdev);
  auto dataNHost = Kokkos::create_mirror_view(dataNdev);
  auto idxHost = Kokkos::create_mirror_view(idxDev);
  auto dataIdxHost = Kokkos::create_mirror_view(dataIdxDev);
  Kokkos::deep_copy(xNHost, xNdev);
  Kokkos::deep_copy(dataNHost, dataNdev);
  Kokkos::deep_copy(idxHost, idxDev);
  Kokkos::deep_copy(dataIdxHost, dataIdxDev);

  //idfx::cout << "2D: x=" << xNHost(0) << "," << xNHost(1) << " y=" << xNHost(2) << ","
  //            << xNHost(3) << " data=" << dataNHost(0) << "," << dataNHost(1) << ","
  //            << dataNHost(2) << "," << dataNHost(3) << " idx=" << idxHost(0) << ","
  //            << idxHost(1) << std::endl;
  real xN[4], dataN[4];
  int idx[2], dataIdx[4];
  for(int n = 0 ; n < 4 ; n++) { xN[n] = xNHost(n); dataN[n] = dataNHost(n); }
  idx[0] = idxHost(0); idx[1] = idxHost(1);
  for(int n = 0 ; n < 4 ; n++) dataIdx[n] = dataIdxHost(n);
  ASSERT_C_ARRAY_EQ(xN, kXN2D, 4, "2D neighbours (device), coordinates");
  ASSERT_C_ARRAY_EQ(dataN, kDataN2D, 4, "2D neighbours (device), data");
  ASSERT_C_ARRAY_EQ(idx, kIdx2D, 2, "2D neighbours (device), idx");
  ASSERT_C_ARRAY_EQ(dataIdx, kDataIdx2D, 4, "2D neighbours (device), dataIdx");
}

//Testing the 1D neighbours on the edges of the table.
TEST(LookupTable, neighbours_on_edges_1d) {
  // Read 1D CSV File
  LookupTable<1> csv1D(DATA_FILE("toto1D.csv"),',');

  // toto1D.csv holds x=1,2,3 and data=2,4,6. Whatever the requested value, we expect the two
  // neighbours bracketing it, i.e. the last two nodes when we sit on the upper edge
  real xq[1];
  real xN[2];
  real dataN[2];
  real expected[3][2] = {{1.0,2.0}, {2.0,3.0}, {2.0,3.0}};
  real xRequest[3] = {1.0, 2.0, 3.0};
  for(int n = 0 ; n < 3 ; n++) {
    xq[0] = xRequest[n];
    csv1D.GetNeighboursHost(xq, xN, dataN);
    //idfx::cout << "x=" << xq[0] << " -> neighbours " << xN[0] << "," << xN[1]
    //            << " (data " << dataN[0] << "," << dataN[1] << ")" << std::endl;
    ASSERT_C_ARRAY_EQ(xN, expected[n], 2, "1D edge neighbours, coordinates");
    EXPECT_LE(xN[0], xq[0]) << "1D edge neighbours do not bracket the requested value";
    EXPECT_GE(xN[1], xq[0]) << "1D edge neighbours do not bracket the requested value";
  }
}

//Testing the reuse of the search between Get and GetNeighbours on Host.
TEST(LookupTable, reuse_search_between_get_and_get_neighbours_on_host) {
  LookupTable<2> csv(DATA_FILE("toto.csv"),',');

  real xq[2];
  xq[0] = 2.1;
  xq[1] = 3.5;
  // The neighbours found by Get are stored in nb...
  LookupTableSearchCache<2> nb;
  real result = csv.GetHost(xq, nb);
  ASSERT_NEAR(result, kValue2D, 1e-13) << "Get-then-GetNeighbours (host), value";
  ASSERT_TRUE(nb.valid) << "Get-then-GetNeighbours (host), search validity";
  ASSERT_C_ARRAY_EQ(nb.idx, kIdx2D, 2, "Get-then-GetNeighbours (host), cached idx");

  // ... and are reused (not computed again) by the getters below
  real xN[4];
  real dataN[4];
  int idx[2];
  int dataIdx[4];
  csv.GetNeighboursHost(xq, nb, xN, dataN);
  csv.GetNeighboursIndxHost(xq, nb, idx, dataIdx);
  ASSERT_C_ARRAY_EQ(xN, kXN2D, 4, "Get-then-GetNeighbours (host), coordinates");
  ASSERT_C_ARRAY_EQ(dataN, kDataN2D, 4, "Get-then-GetNeighbours (host), data");
  ASSERT_C_ARRAY_EQ(idx, kIdx2D, 2, "Get-then-GetNeighbours (host), idx");
  ASSERT_C_ARRAY_EQ(dataIdx, kDataIdx2D, 4, "Get-then-GetNeighbours (host), dataIdx");

  // Check that the table is really not searched again for the same coordinates: we corrupt
  // the stored search, and check that the corrupted result is the one which is used
  nb.idx[0] = 1;
  csv.GetNeighboursIndxHost(xq, nb, idx, dataIdx);
  ASSERT_EQ(idx[0], 1) << "Get-then-GetNeighbours (host): the stored search was not reused";

  // ... while different coordinates trigger a new search, as usual
  real xq2[2];
  xq2[0] = 2.9;
  xq2[1] = 2.5;
  csv.GetNeighboursIndxHost(xq2, nb, idx, dataIdx);
  ASSERT_EQ(idx[0], 0)
    << "Get-then-GetNeighbours (host): different coordinates did not trigger a new search";
  ASSERT_EQ(idx[1], 0)
    << "Get-then-GetNeighbours (host): different coordinates did not trigger a new search";
}

//Testing the search performed by GetNeighbours first, and reused by Get, on Host.
TEST(LookupTable, get_and_reuse_on_host) {
  LookupTable<2> csv(DATA_FILE("toto.csv"),',');
  LookupTable<1> csv1D(DATA_FILE("toto1D.csv"),',');
  real xq[2];
  xq[0] = 2.1;
  xq[1] = 3.5;

  // No call to Get yet: GetNeighbours searches the table as usual, and stores the result
  LookupTableSearchCache<2> nb;
  real xN[4];
  real dataN[4];
  csv.GetNeighboursHost(xq, nb, xN, dataN);
  ASSERT_TRUE(nb.valid) << "GetNeighbours-then-Get (host), search validity";
  ASSERT_C_ARRAY_EQ(nb.idx, kIdx2D, 2, "GetNeighbours-then-Get (host), cached idx");
  ASSERT_C_ARRAY_EQ(xN, kXN2D, 4, "GetNeighbours-then-Get (host), coordinates");
  ASSERT_C_ARRAY_EQ(dataN, kDataN2D, 4, "GetNeighbours-then-Get (host), data");

  // Get now reuses that search for the same coordinates
  real result = csv.GetHost(xq, nb);
  ASSERT_NEAR(result, kValue2D, 1e-13) << "GetNeighbours-then-Get (host), value";

  // Same check as before, the other way around: we corrupt the ratio stored by
  // GetNeighbours, and check that Get uses it instead of computing it again.
  // With delta[0]=0.5 instead of 0.1, the interpolation gives (5+6+6+7)/4=6
  nb.delta[0] = 0.5;
  result = csv.GetHost(xq, nb);
  ASSERT_NEAR(result, 6.0, 1e-13)
    << "GetNeighbours-then-Get (host): the stored search was not reused by Get";

  // ... while different coordinates make Get search the table again
  real xq2[2];
  xq2[0] = 2.9;
  xq2[1] = 2.5;
  result = csv.GetHost(xq2, nb);
  ASSERT_NEAR(result, csv.GetHost(xq2), 1e-13)
    << "GetNeighbours-then-Get (host): different coordinates did not trigger a new search";
  ASSERT_EQ(nb.idx[1], 0)
    << "GetNeighbours-then-Get (host): different coordinates did not trigger a new search";

  // The same holds when the search comes from GetNeighboursIndx, here on the 1D table
  LookupTableSearchCache<1> nb1D;
  real xq1[1];
  xq1[0] = 2.1;
  int idx1[1];
  int dataIdx1[2];
  csv1D.GetNeighboursIndxHost(xq1, nb1D, idx1, dataIdx1);
  ASSERT_TRUE(nb1D.valid) << "GetNeighboursIndx-then-Get (host, 1D), search validity";
  ASSERT_EQ(idx1[0], 1) << "GetNeighboursIndx-then-Get (host, 1D), idx";
  nb1D.delta[0] = 0.0;   // x is now on the left neighbour, so we expect its data
  result = csv1D.GetHost(xq1, nb1D);
  ASSERT_NEAR(result, 4.0, 1e-13)
    << "GetNeighboursIndx-then-Get (host, 1D): the stored search was not reused by Get";
}

//Testing the reuse of the search between Get and GetNeighbours on device.
TEST(LookupTable, get_and_reuse_on_device) {
  IdefixArray1D<real> arr = IdefixArray1D<real>("Test",1);
  IdefixArray1D<real>::host_mirror_type arrHost = Kokkos::create_mirror_view(arr);

  LookupTable<2> csv(DATA_FILE("toto.csv"),',');
  LookupTable<1> csv1D(DATA_FILE("toto1D.csv"),',');

  IdefixArray1D<real> xNdev = IdefixArray1D<real>("xN",4);
  IdefixArray1D<real> dataNdev = IdefixArray1D<real>("dataN",4);
  IdefixArray1D<int> idxDev = IdefixArray1D<int>("idx",2);
  IdefixArray1D<int> dataIdxDev = IdefixArray1D<int>("dataIdx",4);

  idefix_for("neighbours",0, 1, KOKKOS_LAMBDA (int i) {
    real xq[2];
    xq[0] = 2.1;
    xq[1] = 3.5;
    // the structure is local to the loop, and is therefore private to each thread
    LookupTableSearchCache<2> nb;
    real xN[4];
    real dataN[4];
    int idx[2];
    int dataIdx[4];
    arr(i) = csv.Get(xq, nb);
    csv.GetNeighbours(xq, nb, xN, dataN);
    csv.GetNeighboursIndx(xq, nb, idx, dataIdx);
    for(int n = 0 ; n < 4 ; n++) {
      xNdev(n) = xN[n];
      dataNdev(n) = dataN[n];
      dataIdxDev(n) = dataIdx[n];
    }
    idxDev(0) = idx[0];
    idxDev(1) = idx[1];
  });

  Kokkos::deep_copy(arrHost , arr);
  auto xNHost = Kokkos::create_mirror_view(xNdev);
  auto dataNHost = Kokkos::create_mirror_view(dataNdev);
  auto idxHost = Kokkos::create_mirror_view(idxDev);
  auto dataIdxHost = Kokkos::create_mirror_view(dataIdxDev);
  Kokkos::deep_copy(xNHost, xNdev);
  Kokkos::deep_copy(dataNHost, dataNdev);
  Kokkos::deep_copy(idxHost, idxDev);
  Kokkos::deep_copy(dataIdxHost, dataIdxDev);

  // Same test on device, but with GetNeighbours called first and Get reusing its search.
  // The ratio stored by GetNeighbours is corrupted in between, so that a value of 6 (instead
  // of 5.6) proves that Get did not search the table again
  IdefixArray1D<real> valuesDev = IdefixArray1D<real>("values",2);
  idefix_for("neighboursFirst",0, 1, KOKKOS_LAMBDA (int i) {
    real xq[2];
    xq[0] = 2.1;
    xq[1] = 3.5;
    LookupTableSearchCache<2> nb;
    real xN[4];
    real dataN[4];
    // no call to Get yet: the search is performed here for the first time
    csv.GetNeighbours(xq, nb, xN, dataN);
    valuesDev(0) = csv.Get(xq, nb);
    nb.delta[0] = 0.5;
    valuesDev(1) = csv.Get(xq, nb);
  });
  auto valuesHost = Kokkos::create_mirror_view(valuesDev);
  Kokkos::deep_copy(valuesHost, valuesDev);
  //idfx::cout << "neighbours first: value=" << valuesHost(0)
  //            << " (with a corrupted stored search)=" << valuesHost(1) << std::endl;
  ASSERT_NEAR(valuesHost(0), kValue2D, 1e-13)
    << "GetNeighbours-then-Get (device): the stored search was not reused by Get";
  ASSERT_NEAR(valuesHost(1), 6.0, 1e-13)
    << "GetNeighbours-then-Get (device): the stored search was not reused by Get";

  //idfx::cout << "value=" << arrHost(0) << " x=" << xNHost(0) << "," << xNHost(1)
  //            << " y=" << xNHost(2) << "," << xNHost(3)
  //            << " data=" << dataNHost(0) << "," << dataNHost(1) << "," << dataNHost(2)
  //            << "," << dataNHost(3) << " idx=" << idxHost(0) << "," << idxHost(1) << std::endl;
  ASSERT_NEAR(arrHost(0), kValue2D, 1e-13) << "Get-then-GetNeighbours (device), value";
  real xN[4], dataN[4];
  int idx[2], dataIdx[4];
  for(int n = 0 ; n < 4 ; n++) { xN[n] = xNHost(n); dataN[n] = dataNHost(n); }
  idx[0] = idxHost(0); idx[1] = idxHost(1);
  for(int n = 0 ; n < 4 ; n++) dataIdx[n] = dataIdxHost(n);
  ASSERT_C_ARRAY_EQ(xN, kXN2D, 4, "Get-then-GetNeighbours (device), coordinates");
  ASSERT_C_ARRAY_EQ(dataN, kDataN2D, 4, "Get-then-GetNeighbours (device), data");
  ASSERT_C_ARRAY_EQ(idx, kIdx2D, 2, "Get-then-GetNeighbours (device), idx");
  ASSERT_C_ARRAY_EQ(dataIdx, kDataIdx2D, 4, "Get-then-GetNeighbours (device), dataIdx");
}

//Testing 3D npy file on device.
TEST(LookupTable, numpy_file_3d_on_device) {
  IdefixArray1D<real> arr = IdefixArray1D<real>("Test",1);
  IdefixArray1D<real>::host_mirror_type arrHost = Kokkos::create_mirror_view(arr);

  // Read npy File
  std::vector<std::string> coords({
    DATA_FILE("x.npy"),
    DATA_FILE("y.npy"),
    DATA_FILE("z.npy")
  });
  LookupTable<3> csvnpy(coords, std::string(DATA_FILE("data.npy")));

  idefix_for("loop",0, 1, KOKKOS_LAMBDA (int i) {
    real x[3];
    x[0] = 2.7;
    x[1] = 7.4;
    x[2] = 3.9;
    arr(i) = csvnpy.Get(x);
  });

  Kokkos::deep_copy(arrHost , arr);
  //idfx::cout << "result="<<arrHost(0) << std::endl;
  ASSERT_NEAR(arrHost(0), 13.6, 1e-13) << "3D npy, device";
}

//Testing 3D npy file on host.
TEST(LookupTable, numpy_file_3d_on_host) {
  // Read npy File
  std::vector<std::string> coords({
    DATA_FILE("x.npy"),
    DATA_FILE("y.npy"),
    DATA_FILE("z.npy")
  });
  LookupTable<3> csvnpy(coords,
      std::string(DATA_FILE("data.npy"))
  );

  real y[3];
  y[0] = 2.7;
  y[1] = 7.4;
  y[2] = 3.9;
  real result = csvnpy.GetHost(y);

  //idfx::cout << "result="<< result << std::endl;
  ASSERT_NEAR(result, 13.6, 1e-13) << "3D npy, host";
}
