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
    idfx::cout << "--------------------------------------" << std::endl;
    idfx::cout << "Testing 2D CSV file on device." << std::endl;
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
    if(std::fabs(arrHost(0) - 5.6)>1e-13) {
      idfx::cerr << std::scientific;
      idfx::cerr << "ERROR!!" << std::endl;
      idfx::cerr << arrHost(0)-5.6;
      exit(1);
    }
    idfx::cout << "Success" << std::endl;

      idfx::cout << "--------------------------------------" << std::endl;
    idfx::cout << "Testing 2D CSV file on Host." << std::endl;
    real x[2];
    x[0] = 2.1;
    x[1] = 3.5;
    real result = csv.GetHost(x);
    idfx::cout << "result="<<result << std::endl;
    if(std::fabs(result - 5.6)>1e-13) {
      idfx::cerr << std::scientific;
      idfx::cerr << "ERROR!!" << std::endl;
      idfx::cerr << result-5.6;
      exit(1);
    }
    idfx::cout << "Success" << std::endl;

    idfx::cout << "--------------------------------------" << std::endl;
    idfx::cout << "Testing 1D CSV file on device." << std::endl;
    // Read 1D CSV File
    LookupTable<1> csv1D("toto1D.csv",',');

    idefix_for("loop",0, 1, KOKKOS_LAMBDA (int i) {
      real x[2];
      x[0] = 2.1;
      arr(i) = csv1D.Get(x);
    });

    Kokkos::deep_copy(arrHost , arr);

    idfx::cout << "result="<<arrHost(0) << std::endl;
    if(arrHost(0) != 4.2) {
      idfx::cerr << std::scientific;
      idfx::cerr << "ERROR!!" << std::endl;
      exit(1);
    }
    idfx::cout << "Success" << std::endl;

    idfx::cout << "--------------------------------------" << std::endl;
    idfx::cout << "Testing 1D CSV file on Host." << std::endl;

    x[0] = 2.1;
    result = csv1D.GetHost(x);
    idfx::cout << "result="<<result << std::endl;
    if(result != 4.2) {
      idfx::cerr << std::scientific;
      idfx::cerr << "ERROR!!" << std::endl;
      exit(1);
    }
    idfx::cout << "Success" << std::endl;

    idfx::cout << "--------------------------------------" << std::endl;
    idfx::cout << "Testing 1D CSV file read as columns on device." << std::endl;
    // Read the same 1D table, but stored as columns of the CSV file
    LookupTable<1> csv1Dcolumn("toto1Dcolumn.csv",',', true, true);

    idefix_for("loop",0, 1, KOKKOS_LAMBDA (int i) {
      real x[1];
      x[0] = 2.1;
      arr(i) = csv1Dcolumn.Get(x);
    });

    Kokkos::deep_copy(arrHost , arr);

    idfx::cout << "result="<<arrHost(0) << std::endl;
    if(arrHost(0) != 4.2) {
      idfx::cerr << std::scientific;
      idfx::cerr << "ERROR!!" << std::endl;
      exit(1);
    }
    idfx::cout << "Success" << std::endl;

    idfx::cout << "--------------------------------------" << std::endl;
    idfx::cout << "Testing 1D CSV file read as columns on Host." << std::endl;

    x[0] = 2.1;
    result = csv1Dcolumn.GetHost(x);
    idfx::cout << "result="<<result << std::endl;
    if(result != 4.2) {
      idfx::cerr << std::scientific;
      idfx::cerr << "ERROR!!" << std::endl;
      exit(1);
    }
    idfx::cout << "Success" << std::endl;

    idfx::cout << "--------------------------------------" << std::endl;
    idfx::cout << "Testing interpolation in function space on device." << std::endl;
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
    if(std::fabs(arrHost(0) - 9.0)>1e-13) {
      idfx::cerr << std::scientific;
      idfx::cerr << "ERROR!!" << std::endl;
      idfx::cerr << arrHost(0)-9.0;
      exit(1);
    }
    idfx::cout << "Success" << std::endl;

    idefix_for("loop",0, 1, KOKKOS_LAMBDA (int i) {
      real x[1];
      x[0] = 3.0;
      arr(i) = csv1Dlog10.Get(x);
    });

    Kokkos::deep_copy(arrHost , arr);

    idfx::cout << "result (with log10)="<<arrHost(0) << std::endl;
    if(std::fabs(arrHost(0) - 9.0)>1e-13) {
      idfx::cerr << std::scientific;
      idfx::cerr << "ERROR!!" << std::endl;
      idfx::cerr << arrHost(0)-9.0;
      exit(1);
    }
    idfx::cout << "Success" << std::endl;

    idfx::cout << "--------------------------------------" << std::endl;
    idfx::cout << "Testing interpolation in function space on Host." << std::endl;

    x[0] = 3.0;
    result = csv1Dlog.GetHost(x);
    idfx::cout << "result="<<result << std::endl;
    if(std::fabs(result - 9.0)>1e-13) {
      idfx::cerr << std::scientific;
      idfx::cerr << "ERROR!!" << std::endl;
      idfx::cerr << result-9.0;
      exit(1);
    }
    result = csv1Dlog10.GetHost(x);
    idfx::cout << "result (with log10)="<<result << std::endl;
    if(std::fabs(result - 9.0)>1e-13) {
      idfx::cerr << std::scientific;
      idfx::cerr << "ERROR!!" << std::endl;
      idfx::cerr << result-9.0;
      exit(1);
    }
    // the same table, interpolated linearly (default behaviour)
    result = csv1Dlin.GetHost(x);
    idfx::cout << "result (linear interpolation)="<<result << std::endl;
    if(std::fabs(result - 10.0)>1e-13) {
      idfx::cerr << std::scientific;
      idfx::cerr << "ERROR!!" << std::endl;
      idfx::cerr << result-10.0;
      exit(1);
    }
    idfx::cout << "Success" << std::endl;

    idfx::cout << "--------------------------------------" << std::endl;
    idfx::cout << "Testing the neighbours used for the interpolation on Host." << std::endl;
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
      if(xN[0] != 2.0 || xN[1] != 3.0 || dataN[0] != 4.0 || dataN[1] != 6.0
         || idx[0] != 1 || dataIdx[0] != 1 || dataIdx[1] != 2) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }

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
      if(xN2[0] != 2.0 || xN2[1] != 3.0 || xN2[2] != 3.0 || xN2[3] != 4.0) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }
      if(dataN2[0] != 5.0 || dataN2[1] != 6.0 || dataN2[2] != 6.0 || dataN2[3] != 7.0) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }
      if(idx2[0] != 0 || idx2[1] != 1
         || dataIdx2[0] != 1 || dataIdx2[1] != 4 || dataIdx2[2] != 2 || dataIdx2[3] != 5) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }

      // The neighbours should reproduce the value returned by GetHost
      real dx = (xq2[0]-xN2[0])/(xN2[1]-xN2[0]);
      real dy = (xq2[1]-xN2[2])/(xN2[3]-xN2[2]);
      real interpolated = (1-dx)*(1-dy)*dataN2[0] + dx*(1-dy)*dataN2[1]
                        + (1-dx)*dy*dataN2[2] + dx*dy*dataN2[3];
      idfx::cout << "2D: interpolation from the neighbours=" << interpolated << std::endl;
      if(std::fabs(interpolated - csv.GetHost(xq2))>1e-13) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }

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
      if(std::fabs(xN3[0]-2.0)>1e-13 || std::fabs(xN3[1]-4.0)>1e-13
         || std::fabs(dataN3[0]-4.0)>1e-13 || std::fabs(dataN3[1]-16.0)>1e-13
         || idx3[0] != 1 || dataIdx3[0] != 1 || dataIdx3[1] != 2) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }
      idfx::cout << "Success" << std::endl;
    }

    idfx::cout << "--------------------------------------" << std::endl;
    idfx::cout << "Testing the neighbours used for the interpolation on device." << std::endl;
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
      if(xNHost(0) != 2.0 || xNHost(1) != 3.0 || xNHost(2) != 3.0 || xNHost(3) != 4.0) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }
      if(dataNHost(0) != 5.0 || dataNHost(1) != 6.0 || dataNHost(2) != 6.0
         || dataNHost(3) != 7.0) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }
      if(idxHost(0) != 0 || idxHost(1) != 1 || dataIdxHost(0) != 1 || dataIdxHost(1) != 4
         || dataIdxHost(2) != 2 || dataIdxHost(3) != 5) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }
      idfx::cout << "Success" << std::endl;
    }

    idfx::cout << "--------------------------------------" << std::endl;
    idfx::cout << "Testing the 1D neighbours on the edges of the table." << std::endl;
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
        if(xN[0] != expected[n][0] || xN[1] != expected[n][1]) {
          idfx::cerr << "ERROR!!" << std::endl;
          exit(1);
        }
        if(xN[0] > xq[0] || xN[1] < xq[0]) {
          idfx::cerr << "ERROR!! the neighbours do not bracket the requested value" << std::endl;
          exit(1);
        }
      }
      idfx::cout << "Success" << std::endl;
    }

    idfx::cout << "--------------------------------------" << std::endl;
    idfx::cout << "Testing the reuse of the search between Get and GetNeighbours on Host."
               << std::endl;
    {
      real xq[2];
      xq[0] = 2.1;
      xq[1] = 3.5;
      // The neighbours found by Get are stored in nb...
      LookupTableNeighbours<2> nb;
      result = csv.GetHost(xq, nb);
      if(std::fabs(result - 5.6)>1e-13 || !nb.valid || nb.idx[0] != 0 || nb.idx[1] != 1) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }

      // ... and are reused (not computed again) by the getters below
      real xN[4];
      real dataN[4];
      int idx[2];
      int dataIdx[4];
      csv.GetNeighboursHost(xq, nb, xN, dataN);
      csv.GetNeighboursIndxHost(xq, nb, idx, dataIdx);
      if(xN[0] != 2.0 || xN[1] != 3.0 || xN[2] != 3.0 || xN[3] != 4.0) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }
      if(dataN[0] != 5.0 || dataN[1] != 6.0 || dataN[2] != 6.0 || dataN[3] != 7.0) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }
      if(idx[0] != 0 || idx[1] != 1
         || dataIdx[0] != 1 || dataIdx[1] != 4 || dataIdx[2] != 2 || dataIdx[3] != 5) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }

      // Check that the table is really not searched again for the same coordinates: we corrupt
      // the stored search, and check that the corrupted result is the one which is used
      nb.idx[0] = 1;
      csv.GetNeighboursIndxHost(xq, nb, idx, dataIdx);
      if(idx[0] != 1) {
        idfx::cerr << "ERROR!! the stored search was not reused" << std::endl;
        exit(1);
      }

      // ... while different coordinates trigger a new search, as usual
      real xq2[2];
      xq2[0] = 2.9;
      xq2[1] = 2.5;
      csv.GetNeighboursIndxHost(xq2, nb, idx, dataIdx);
      if(idx[0] != 0 || idx[1] != 0) {
        idfx::cerr << "ERROR!! the search was not performed again" << std::endl;
        exit(1);
      }
      idfx::cout << "Success" << std::endl;
    }

    idfx::cout << "--------------------------------------" << std::endl;
    idfx::cout << "Testing the search performed by GetNeighbours first, and reused by Get, "
                  "on Host." << std::endl;
    {
      real xq[2];
      xq[0] = 2.1;
      xq[1] = 3.5;

      // No call to Get yet: GetNeighbours searches the table as usual, and stores the result
      LookupTableNeighbours<2> nb;
      real xN[4];
      real dataN[4];
      csv.GetNeighboursHost(xq, nb, xN, dataN);
      if(!nb.valid || nb.idx[0] != 0 || nb.idx[1] != 1) {
        idfx::cerr << "ERROR!! the search was not performed" << std::endl;
        exit(1);
      }
      if(xN[0] != 2.0 || xN[1] != 3.0 || dataN[0] != 5.0 || dataN[3] != 7.0) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }

      // Get now reuses that search for the same coordinates
      result = csv.GetHost(xq, nb);
      if(std::fabs(result - 5.6)>1e-13) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }

      // Same check as before, the other way around: we corrupt the ratio stored by
      // GetNeighbours, and check that Get uses it instead of computing it again.
      // With delta[0]=0.5 instead of 0.1, the interpolation gives (5+6+6+7)/4=6
      nb.delta[0] = 0.5;
      result = csv.GetHost(xq, nb);
      if(std::fabs(result - 6.0)>1e-13) {
        idfx::cerr << "ERROR!! the stored search was not reused by Get" << std::endl;
        exit(1);
      }

      // ... while different coordinates make Get search the table again
      real xq2[2];
      xq2[0] = 2.9;
      xq2[1] = 2.5;
      result = csv.GetHost(xq2, nb);
      if(std::fabs(result - csv.GetHost(xq2))>1e-13 || nb.idx[1] != 0) {
        idfx::cerr << "ERROR!! the search was not performed again" << std::endl;
        exit(1);
      }

      // The same holds when the search comes from GetNeighboursIndx, here on the 1D table
      LookupTableNeighbours<1> nb1D;
      real xq1[1];
      xq1[0] = 2.1;
      int idx1[1];
      int dataIdx1[2];
      csv1D.GetNeighboursIndxHost(xq1, nb1D, idx1, dataIdx1);
      if(!nb1D.valid || idx1[0] != 1) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }
      nb1D.delta[0] = 0.0;   // x is now on the left neighbour, so we expect its data
      result = csv1D.GetHost(xq1, nb1D);
      if(std::fabs(result - 4.0)>1e-13) {
        idfx::cerr << "ERROR!! the stored search was not reused by Get" << std::endl;
        exit(1);
      }
      idfx::cout << "Success" << std::endl;
    }

    idfx::cout << "--------------------------------------" << std::endl;
    idfx::cout << "Testing the reuse of the search between Get and GetNeighbours on device."
               << std::endl;
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
        LookupTableNeighbours<2> nb;
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
        LookupTableNeighbours<2> nb;
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
      if(std::fabs(valuesHost(0) - 5.6)>1e-13 || std::fabs(valuesHost(1) - 6.0)>1e-13) {
        idfx::cerr << "ERROR!! the stored search was not reused by Get" << std::endl;
        exit(1);
      }

      idfx::cout << "value=" << arrHost(0) << " x=" << xNHost(0) << "," << xNHost(1)
                 << " y=" << xNHost(2) << "," << xNHost(3)
                 << " data=" << dataNHost(0) << "," << dataNHost(1) << "," << dataNHost(2)
                 << "," << dataNHost(3) << " idx=" << idxHost(0) << "," << idxHost(1) << std::endl;
      if(std::fabs(arrHost(0) - 5.6)>1e-13) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }
      if(xNHost(0) != 2.0 || xNHost(1) != 3.0 || xNHost(2) != 3.0 || xNHost(3) != 4.0) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }
      if(dataNHost(0) != 5.0 || dataNHost(1) != 6.0 || dataNHost(2) != 6.0
         || dataNHost(3) != 7.0) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }
      if(idxHost(0) != 0 || idxHost(1) != 1 || dataIdxHost(0) != 1 || dataIdxHost(1) != 4
         || dataIdxHost(2) != 2 || dataIdxHost(3) != 5) {
        idfx::cerr << "ERROR!!" << std::endl;
        exit(1);
      }
      idfx::cout << "Success" << std::endl;
    }

    idfx::cout << "--------------------------------------" << std::endl;
    idfx::cout << "Testing 3D npy file on device." << std::endl;
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
    if(std::fabs(arrHost(0) - 13.6)>1e-13) {
      idfx::cerr << std::scientific;
      idfx::cerr << "ERROR!!" << std::endl;
      idfx::cerr << arrHost(0)-13.6;
      exit(1);
    }
    idfx::cout << "Success" << std::endl;

    idfx::cout << "--------------------------------------" << std::endl;
    idfx::cout << "Testing 3D npy file on host." << std::endl;

    real y[3];
    y[0] = 2.7;
    y[1] = 7.4;
    y[2] = 3.9;
    result = csvnpy.GetHost(y);

    idfx::cout << "result="<< result << std::endl;
    if(std::fabs(result - 13.6)>1e-13) {
      idfx::cerr << std::scientific;
      idfx::cerr << "ERROR!!" << std::endl;
      idfx::cerr << result-13.6;
      exit(1);
    }
    idfx::cout << "Success" << std::endl;
    idfx::cout << "--------------------------------------" << std::endl;
    idfx::cout << "Done." << std::endl;

  }
  Kokkos::finalize();
  #ifdef WITH_MPI
    MPI_Finalize();
  #endif

  return 0;

}
