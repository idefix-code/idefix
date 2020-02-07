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







int main( int argc, char* argv[] )
{

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
