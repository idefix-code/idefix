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
    GridHost gridHost(grid);

    // Actually make the grid
    gridHost.MakeGrid(input);
    gridHost.SyncToDevice();

    // Make a data array
    DataBlock data;

    data.InitFromGrid(grid);

    DataBlockHost dataHost(data);
    dataHost.SyncFromDevice();

    std::cout << "Init Physics" << std::endl;
    Physics phys(input);
    phys.InitFlow(data);

    std::cout << "Init Time Integrator..." << std::endl;
    TimeIntegrator Tint(input, phys);

    std::cout << "init Output Routines" << std::endl;
    OutputVTK output = OutputVTK(input, grid, Tint.getT()); 

    std::cout << "Write init vtk" << std::endl;
    output.Write(data,Tint.getT());
    std::cout << "Cycling Time Integrator..." << std::endl;

    Kokkos::Timer timer;

    while(Tint.getT() < input.tfinal) {
      Tint.Cycle(data);
      output.Write(data, Tint.getT());
    }
    double tintegration = (timer.seconds()/(grid.np_int[IDIR]*grid.np_int[JDIR]*grid.np_int[KDIR]*Tint.getNcycles()));

    std::cout << "Completed in " << timer.seconds() << "seconds. Perfs are " << 1/tintegration << " cell updates/second." << std::endl;
    

    // Make a test
    /*
    Test test(data);
    int nrepeat=10000;
    test.MakeTest(IDIR-1,nrepeat);
    test.MakeTest(IDIR,nrepeat);
    test.MakeTest(JDIR,nrepeat);
    test.MakeTest(KDIR,nrepeat);
    */
    std::cout << "Job's done" << std::endl;
  }
  Kokkos::finalize();

  return 0;
}
