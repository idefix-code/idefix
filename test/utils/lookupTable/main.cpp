#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <Kokkos_Core.hpp>

#include "idefix.hpp"
#include "lookupTable.hpp"

// minimal skeleton to use idfx basic functions
void testReduction();


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
    IdefixArray1D<real>::HostMirror arrHost = Kokkos::create_mirror_view(arr);

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
