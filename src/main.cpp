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
    DataBlock data;

    data.InitFromGrid(grid);
    DataBlockHost dataHost(data);

    dataHost.SyncFromDevice();

    std::cout << "Init Time Integrator..." << std::endl;
    TimeIntegrator Tint(input, data);
    Tint.setDt(0.1);

    std::cout << "Cycling Time Integrator..." << std::endl;
    Tint.Cycle();

    std::cout << "Done." << std::endl;

    // Make a test
    Test test(data);
    int nrepeat=10000;
    test.MakeTest(IDIR-1,nrepeat);
    test.MakeTest(IDIR,nrepeat);
    test.MakeTest(JDIR,nrepeat);
    test.MakeTest(KDIR,nrepeat);

    std::cout << "Job's done" << std::endl;
  }
  Kokkos::finalize();

  return 0;
}
