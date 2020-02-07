/*
//@HEADER
// ************************************************************************
//
//                        IDEFIX v 0.0-alpha
//
// ************************************************************************
//@HEADER
*/

#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>



#include <Kokkos_Core.hpp>

#include "idefix.hpp"

#define  I_STRIDE
//#define  J_STRIDE
//#define  K_STRIDE
//#define  NO_STRIDE



//#define USE_LEFT_ITERATE





int main( int argc, char* argv[] )
{
  const int NX=64;
  const int NY=64;
  const int NZ=64;

#ifdef I_STRIDE
  const int is=1;
  const int ie=NX-2;
#else
  const int is=0;
  const int ie=NX-1;
#endif

#ifdef J_STRIDE
  const int js=1;
  const int je=NY-2;
#else
  const int js=0;
  const int je=NY-1;
#endif

#ifdef K_STRIDE
  const int ks=1;
  const int ke=NZ-2;
#else
  const int ks=0;
  const int ke=NZ-1;
#endif


  int nrepeat=10000;

  Kokkos::initialize( argc, argv );
  {

    Input input = Input();

    // Allocate the grid on device
    Grid grid(input);

    // Allocate the grid on host as a mirror
    GridHost gridHost(grid,input);

    // Actually make the grid
    gridHost.MakeGrid(grid,input);

    // Make a data array
    Data data(grid);
    DataHost dataHost(data);

    dataHost.SyncFromDevice();

    // Make a test
    Test test(data);
    int nrepeat=10000;
    test.MakeTest(grid,IDIR-1,nrepeat);
    test.MakeTest(grid,IDIR,nrepeat);
    test.MakeTest(grid,JDIR,nrepeat);
    test.MakeTest(grid,KDIR,nrepeat);

  }
  Kokkos::finalize();

  return 0;
}
