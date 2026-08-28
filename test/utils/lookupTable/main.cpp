#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <Kokkos_Core.hpp>

#include "idefix.hpp"
#include "lookupTable.hpp"

// minimal skeleton to use idfx basic functions
void testReduction();

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

// ---------------------------------------------------------------------------------------------
// Small test helpers, so that each check below states what it verifies and why it failed, instead
// of every failure printing an undifferentiated "ERROR!!" that has to be traced back to the last
// banner that was printed.
// ---------------------------------------------------------------------------------------------
namespace {

void Banner(const std::string &what) {
  idfx::cout << "--------------------------------------" << std::endl;
  idfx::cout << what << std::endl;
}

void Success() {
  idfx::cout << "Success" << std::endl;
}

// Exact equality: used for values that should reproduce bit-for-bit (e.g. requesting a table node
// directly, or comparing an int index), matching what the original tests checked with a plain !=.
template<typename T>
void CheckEqual(T got, T expected, const std::string &what) {
  if(got != expected) {
    idfx::cerr << "ERROR!! " << what << ": got " << got << ", expected " << expected << std::endl;
    exit(1);
  }
}

// Element-wise exact equality over a small fixed-size array (neighbour coordinates/data/indices).
template<typename T>
void CheckArrayEqual(const T *got, const T *expected, int n, const std::string &what) {
  for(int i = 0 ; i < n ; i++) {
    if(got[i] != expected[i]) {
      idfx::cerr << "ERROR!! " << what << " (element " << i << "): got " << got[i]
                 << ", expected " << expected[i] << std::endl;
      exit(1);
    }
  }
}

// Tolerance-based equality, for values obtained through interpolation.
void CheckClose(real got, real expected, real tol, const std::string &what) {
  if(std::fabs(got - expected) > tol) {
    idfx::cerr << std::scientific;
    idfx::cerr << "ERROR!! " << what << ": got " << got << ", expected " << expected
               << " (|diff|=" << std::fabs(got - expected) << ", tol=" << tol << ")" << std::endl;
    exit(1);
  }
}

void CheckTrue(bool cond, const std::string &what) {
  if(!cond) {
    idfx::cerr << "ERROR!! " << what << std::endl;
    exit(1);
  }
}

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

}  // namespace

// main function
int main( int argc, char* argv[] )
{
  bool initKokkosBeforeMPI = false;

  // When running on GPUS with Omnipath network,
  // Kokkos needs to be initialised *before* the MPI layer
#ifdef KOKKOS_ENABLE_CUDA
  if(std::getenv("PSM2_CUDA") != NULL) {
    initKokkosBeforeMPI = true;
  }
#endif

  if(initKokkosBeforeMPI)  Kokkos::initialize( argc, argv );

#ifdef WITH_MPI
  MPI_Init(&argc,&argv);
#endif

  if(!initKokkosBeforeMPI) Kokkos::initialize( argc, argv );


  {
    idfx::initialize();

    Banner("Testing 2D CSV file on device.");
    IdefixArray1D<real> arr = IdefixArray1D<real>("Test",1);
    IdefixArray1D<real>::host_mirror_type arrHost = Kokkos::create_mirror_view(arr);

    LookupTable<2> csv("toto.csv",',');

    idefix_for("loop",0, 1, KOKKOS_LAMBDA (int i) {
      real x[2];
      x[0] = 2.1;
      x[1] = 3.5;
      arr(i) = csv.Get(x);
    });

    Kokkos::deep_copy(arrHost , arr);
    idfx::cout << "result="<<arrHost(0) << std::endl;
    CheckClose(arrHost(0), kValue2D, 1e-13, "2D CSV, device");
    Success();

    Banner("Testing 2D CSV file on Host.");
    real x[2];
    x[0] = 2.1;
    x[1] = 3.5;
    real result = csv.GetHost(x);
    idfx::cout << "result="<<result << std::endl;
    CheckClose(result, kValue2D, 1e-13, "2D CSV, host");
    Success();

    Banner("Testing 1D CSV file on device.");
    // Read 1D CSV File
    LookupTable<1> csv1D("toto1D.csv",',');

    idefix_for("loop",0, 1, KOKKOS_LAMBDA (int i) {
      real x[2];
      x[0] = 2.1;
      arr(i) = csv1D.Get(x);
    });

    Kokkos::deep_copy(arrHost , arr);
    idfx::cout << "result="<<arrHost(0) << std::endl;
    CheckEqual(arrHost(0), real(4.2), "1D CSV, device");
    Success();

    Banner("Testing 1D CSV file on Host.");
    x[0] = 2.1;
    result = csv1D.GetHost(x);
    idfx::cout << "result="<<result << std::endl;
    CheckEqual(result, real(4.2), "1D CSV, host");
    Success();

    Banner("Testing 1D CSV file read as columns on device.");
    // Read the same 1D table, but stored as columns of the CSV file
    LookupTable<1> csv1Dcolumn("toto1Dcolumn.csv",',', true, true);

    idefix_for("loop",0, 1, KOKKOS_LAMBDA (int i) {
      real x[1];
      x[0] = 2.1;
      arr(i) = csv1Dcolumn.Get(x);
    });

    Kokkos::deep_copy(arrHost , arr);
    idfx::cout << "result="<<arrHost(0) << std::endl;
    CheckEqual(arrHost(0), real(4.2), "1D CSV read as columns, device");
    Success();

    Banner("Testing 1D CSV file read as columns on Host.");
    x[0] = 2.1;
    result = csv1Dcolumn.GetHost(x);
    idfx::cout << "result="<<result << std::endl;
    CheckEqual(result, real(4.2), "1D CSV read as columns, host");
    Success();

    Banner("Testing interpolation in function space on device.");
    // data = x^2 sampled in x=1,2,4. A linear interpolation in x=3 gives 10, while an
    // interpolation performed on log(x) and log(data) gives the exact result 9
    LookupTable<1> csv1Dlin("toto1Dlog.csv",',', true, true);
    LookupTable<1> csv1Dlog("toto1Dlog.csv",',', true, true, true);
    // the transformation and its inverse can also be chosen when the table is created
    LookupTable<1, MyLog10, MyPow10> csv1Dlog10("toto1Dlog.csv",',', true, true, true);

    idefix_for("loop",0, 1, KOKKOS_LAMBDA (int i) {
      real x[1];
      x[0] = 3.0;
      arr(i) = csv1Dlog.Get(x);
    });
    Kokkos::deep_copy(arrHost , arr);
    idfx::cout << "result="<<arrHost(0) << std::endl;
    CheckClose(arrHost(0), 9.0, 1e-13, "function-space interpolation (default log), device");

    idefix_for("loop",0, 1, KOKKOS_LAMBDA (int i) {
      real x[1];
      x[0] = 3.0;
      arr(i) = csv1Dlog10.Get(x);
    });
    Kokkos::deep_copy(arrHost , arr);
    idfx::cout << "result (with log10)="<<arrHost(0) << std::endl;
    CheckClose(arrHost(0), 9.0, 1e-13, "function-space interpolation (custom log10), device");
    Success();

    Banner("Testing interpolation in function space on Host.");
    x[0] = 3.0;
    result = csv1Dlog.GetHost(x);
    idfx::cout << "result="<<result << std::endl;
    CheckClose(result, 9.0, 1e-13, "function-space interpolation (default log), host");

    result = csv1Dlog10.GetHost(x);
    idfx::cout << "result (with log10)="<<result << std::endl;
    CheckClose(result, 9.0, 1e-13, "function-space interpolation (custom log10), host");

    // the same table, interpolated linearly (default behaviour)
    result = csv1Dlin.GetHost(x);
    idfx::cout << "result (linear interpolation)="<<result << std::endl;
    CheckClose(result, 10.0, 1e-13, "plain linear interpolation, host");
    Success();

    Banner("Testing the neighbours used for the interpolation on Host.");
    {
      // 1D table (x=1,2,3 and data=2,4,6) interpolated in x=2.1
      real xq[1];
      xq[0] = 2.1;
      real xN[2];
      real dataN[2];
      int idx[1];
      int dataIdx[2];
      csv1D.GetNeighboursHost(xq, xN, dataN);
      csv1D.GetNeighboursIndxHost(xq, idx, dataIdx);
      idfx::cout << "1D: x=" << xN[0] << "," << xN[1] << " data=" << dataN[0] << "," << dataN[1]
                 << " idx=" << idx[0] << " dataIdx=" << dataIdx[0] << "," << dataIdx[1]
                 << std::endl;
      const real expectedXN[2] = {2.0, 3.0};
      const real expectedDataN[2] = {4.0, 6.0};
      const int expectedDataIdx[2] = {1, 2};
      CheckArrayEqual(xN, expectedXN, 2, "1D neighbours, coordinates");
      CheckArrayEqual(dataN, expectedDataN, 2, "1D neighbours, data");
      CheckEqual(idx[0], 1, "1D neighbours, idx");
      CheckArrayEqual(dataIdx, expectedDataIdx, 2, "1D neighbours, dataIdx");

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
      idfx::cout << "2D: x=" << xN2[0] << "," << xN2[1] << " y=" << xN2[2] << "," << xN2[3]
                 << " data=" << dataN2[0] << "," << dataN2[1] << "," << dataN2[2] << ","
                 << dataN2[3] << " idx=" << idx2[0] << "," << idx2[1] << std::endl;
      CheckArrayEqual(xN2, kXN2D, 4, "2D neighbours, coordinates");
      CheckArrayEqual(dataN2, kDataN2D, 4, "2D neighbours, data");
      CheckArrayEqual(idx2, kIdx2D, 2, "2D neighbours, idx");
      CheckArrayEqual(dataIdx2, kDataIdx2D, 4, "2D neighbours, dataIdx");

      // The neighbours should reproduce the value returned by GetHost
      real dx = (xq2[0]-xN2[0])/(xN2[1]-xN2[0]);
      real dy = (xq2[1]-xN2[2])/(xN2[3]-xN2[2]);
      real interpolated = (1-dx)*(1-dy)*dataN2[0] + dx*(1-dy)*dataN2[1]
                        + (1-dx)*dy*dataN2[2] + dx*dy*dataN2[3];
      idfx::cout << "2D: interpolation from the neighbours=" << interpolated << std::endl;
      CheckClose(interpolated, csv.GetHost(xq2), 1e-13,
                 "2D neighbours, manual interpolation vs GetHost");

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
      idfx::cout << "1D (function space): x=" << xN3[0] << "," << xN3[1]
                 << " data=" << dataN3[0] << "," << dataN3[1] << " idx=" << idx3[0] << std::endl;
      CheckClose(xN3[0], 2.0, 1e-13, "1D function-space neighbours, x lower bound");
      CheckClose(xN3[1], 4.0, 1e-13, "1D function-space neighbours, x upper bound");
      CheckClose(dataN3[0], 4.0, 1e-13, "1D function-space neighbours, data lower bound");
      CheckClose(dataN3[1], 16.0, 1e-13, "1D function-space neighbours, data upper bound");
      CheckEqual(idx3[0], 1, "1D function-space neighbours, idx");
      const int expectedDataIdx3[2] = {1, 2};
      CheckArrayEqual(dataIdx3, expectedDataIdx3, 2, "1D function-space neighbours, dataIdx");
      Success();
    }

    Banner("Testing the neighbours used for the interpolation on device.");
    {
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

      idfx::cout << "2D: x=" << xNHost(0) << "," << xNHost(1) << " y=" << xNHost(2) << ","
                 << xNHost(3) << " data=" << dataNHost(0) << "," << dataNHost(1) << ","
                 << dataNHost(2) << "," << dataNHost(3) << " idx=" << idxHost(0) << ","
                 << idxHost(1) << std::endl;
      real xN[4], dataN[4];
      int idx[2], dataIdx[4];
      for(int n = 0 ; n < 4 ; n++) { xN[n] = xNHost(n); dataN[n] = dataNHost(n); }
      idx[0] = idxHost(0); idx[1] = idxHost(1);
      for(int n = 0 ; n < 4 ; n++) dataIdx[n] = dataIdxHost(n);
      CheckArrayEqual(xN, kXN2D, 4, "2D neighbours (device), coordinates");
      CheckArrayEqual(dataN, kDataN2D, 4, "2D neighbours (device), data");
      CheckArrayEqual(idx, kIdx2D, 2, "2D neighbours (device), idx");
      CheckArrayEqual(dataIdx, kDataIdx2D, 4, "2D neighbours (device), dataIdx");
      Success();
    }

    Banner("Testing the 1D neighbours on the edges of the table.");
    {
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
        idfx::cout << "x=" << xq[0] << " -> neighbours " << xN[0] << "," << xN[1]
                   << " (data " << dataN[0] << "," << dataN[1] << ")" << std::endl;
        CheckArrayEqual(xN, expected[n], 2, "1D edge neighbours, coordinates");
        CheckTrue(xN[0] <= xq[0] && xN[1] >= xq[0],
                  "1D edge neighbours do not bracket the requested value");
      }
      Success();
    }

    Banner("Testing the reuse of the search between Get and GetNeighbours on Host.");
    {
      real xq[2];
      xq[0] = 2.1;
      xq[1] = 3.5;
      // The neighbours found by Get are stored in nb...
      LookupTableSearchCache<2> nb;
      result = csv.GetHost(xq, nb);
      CheckClose(result, kValue2D, 1e-13, "Get-then-GetNeighbours (host), value");
      CheckTrue(nb.valid, "Get-then-GetNeighbours (host), search validity");
      CheckArrayEqual(nb.idx, kIdx2D, 2, "Get-then-GetNeighbours (host), cached idx");

      // ... and are reused (not computed again) by the getters below
      real xN[4];
      real dataN[4];
      int idx[2];
      int dataIdx[4];
      csv.GetNeighboursHost(xq, nb, xN, dataN);
      csv.GetNeighboursIndxHost(xq, nb, idx, dataIdx);
      CheckArrayEqual(xN, kXN2D, 4, "Get-then-GetNeighbours (host), coordinates");
      CheckArrayEqual(dataN, kDataN2D, 4, "Get-then-GetNeighbours (host), data");
      CheckArrayEqual(idx, kIdx2D, 2, "Get-then-GetNeighbours (host), idx");
      CheckArrayEqual(dataIdx, kDataIdx2D, 4, "Get-then-GetNeighbours (host), dataIdx");

      // Check that the table is really not searched again for the same coordinates: we corrupt
      // the stored search, and check that the corrupted result is the one which is used
      nb.idx[0] = 1;
      csv.GetNeighboursIndxHost(xq, nb, idx, dataIdx);
      CheckEqual(idx[0], 1, "Get-then-GetNeighbours (host): the stored search was not reused");

      // ... while different coordinates trigger a new search, as usual
      real xq2[2];
      xq2[0] = 2.9;
      xq2[1] = 2.5;
      csv.GetNeighboursIndxHost(xq2, nb, idx, dataIdx);
      CheckTrue(idx[0] == 0 && idx[1] == 0,
                "Get-then-GetNeighbours (host): different coordinates did not trigger a new search");
      Success();
    }

    Banner("Testing the search performed by GetNeighbours first, and reused by Get, on Host.");
    {
      real xq[2];
      xq[0] = 2.1;
      xq[1] = 3.5;

      // No call to Get yet: GetNeighbours searches the table as usual, and stores the result
      LookupTableSearchCache<2> nb;
      real xN[4];
      real dataN[4];
      csv.GetNeighboursHost(xq, nb, xN, dataN);
      CheckTrue(nb.valid, "GetNeighbours-then-Get (host), search validity");
      CheckArrayEqual(nb.idx, kIdx2D, 2, "GetNeighbours-then-Get (host), cached idx");
      CheckArrayEqual(xN, kXN2D, 4, "GetNeighbours-then-Get (host), coordinates");
      CheckArrayEqual(dataN, kDataN2D, 4, "GetNeighbours-then-Get (host), data");

      // Get now reuses that search for the same coordinates
      result = csv.GetHost(xq, nb);
      CheckClose(result, kValue2D, 1e-13, "GetNeighbours-then-Get (host), value");

      // Same check as before, the other way around: we corrupt the ratio stored by
      // GetNeighbours, and check that Get uses it instead of computing it again.
      // With delta[0]=0.5 instead of 0.1, the interpolation gives (5+6+6+7)/4=6
      nb.delta[0] = 0.5;
      result = csv.GetHost(xq, nb);
      CheckClose(result, 6.0, 1e-13,
                 "GetNeighbours-then-Get (host): the stored search was not reused by Get");

      // ... while different coordinates make Get search the table again
      real xq2[2];
      xq2[0] = 2.9;
      xq2[1] = 2.5;
      result = csv.GetHost(xq2, nb);
      CheckClose(result, csv.GetHost(xq2), 1e-13,
                 "GetNeighbours-then-Get (host): different coordinates did not trigger a new search");
      CheckEqual(nb.idx[1], 0,
                 "GetNeighbours-then-Get (host): different coordinates did not trigger a new search");

      // The same holds when the search comes from GetNeighboursIndx, here on the 1D table
      LookupTableSearchCache<1> nb1D;
      real xq1[1];
      xq1[0] = 2.1;
      int idx1[1];
      int dataIdx1[2];
      csv1D.GetNeighboursIndxHost(xq1, nb1D, idx1, dataIdx1);
      CheckTrue(nb1D.valid, "GetNeighboursIndx-then-Get (host, 1D), search validity");
      CheckEqual(idx1[0], 1, "GetNeighboursIndx-then-Get (host, 1D), idx");
      nb1D.delta[0] = 0.0;   // x is now on the left neighbour, so we expect its data
      result = csv1D.GetHost(xq1, nb1D);
      CheckClose(result, 4.0, 1e-13,
                 "GetNeighboursIndx-then-Get (host, 1D): the stored search was not reused by Get");
      Success();
    }

    Banner("Testing the reuse of the search between Get and GetNeighbours on device.");
    {
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
      idfx::cout << "neighbours first: value=" << valuesHost(0)
                 << " (with a corrupted stored search)=" << valuesHost(1) << std::endl;
      CheckClose(valuesHost(0), kValue2D, 1e-13,
                 "GetNeighbours-then-Get (device): the stored search was not reused by Get");
      CheckClose(valuesHost(1), 6.0, 1e-13,
                 "GetNeighbours-then-Get (device): the stored search was not reused by Get");

      idfx::cout << "value=" << arrHost(0) << " x=" << xNHost(0) << "," << xNHost(1)
                 << " y=" << xNHost(2) << "," << xNHost(3)
                 << " data=" << dataNHost(0) << "," << dataNHost(1) << "," << dataNHost(2)
                 << "," << dataNHost(3) << " idx=" << idxHost(0) << "," << idxHost(1) << std::endl;
      CheckClose(arrHost(0), kValue2D, 1e-13, "Get-then-GetNeighbours (device), value");
      real xN[4], dataN[4];
      int idx[2], dataIdx[4];
      for(int n = 0 ; n < 4 ; n++) { xN[n] = xNHost(n); dataN[n] = dataNHost(n); }
      idx[0] = idxHost(0); idx[1] = idxHost(1);
      for(int n = 0 ; n < 4 ; n++) dataIdx[n] = dataIdxHost(n);
      CheckArrayEqual(xN, kXN2D, 4, "Get-then-GetNeighbours (device), coordinates");
      CheckArrayEqual(dataN, kDataN2D, 4, "Get-then-GetNeighbours (device), data");
      CheckArrayEqual(idx, kIdx2D, 2, "Get-then-GetNeighbours (device), idx");
      CheckArrayEqual(dataIdx, kDataIdx2D, 4, "Get-then-GetNeighbours (device), dataIdx");
      Success();
    }

    Banner("Testing 3D npy file on device.");
    // Read npy File
    std::vector<std::string> coords({"x.npy","y.npy","z.npy"});

    LookupTable<3> csvnpy(coords, std::string("data.npy"));

    idefix_for("loop",0, 1, KOKKOS_LAMBDA (int i) {
      real x[3];
      x[0] = 2.7;
      x[1] = 7.4;
      x[2] = 3.9;
      arr(i) = csvnpy.Get(x);
    });

    Kokkos::deep_copy(arrHost , arr);
    idfx::cout << "result="<<arrHost(0) << std::endl;
    CheckClose(arrHost(0), 13.6, 1e-13, "3D npy, device");
    Success();

    Banner("Testing 3D npy file on host.");
    real y[3];
    y[0] = 2.7;
    y[1] = 7.4;
    y[2] = 3.9;
    result = csvnpy.GetHost(y);
    idfx::cout << "result="<< result << std::endl;
    CheckClose(result, 13.6, 1e-13, "3D npy, host");
    Success();

    Banner("Done.");
  }
  Kokkos::finalize();
  #ifdef WITH_MPI
    MPI_Finalize();
  #endif

  return 0;

}
