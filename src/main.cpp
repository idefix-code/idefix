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

#ifdef WITH_MPI    
  MPI_Init(&argc,&argv);
#endif

  Kokkos::initialize( argc, argv );
  {

    idfx::initialize();

    Input input = Input("idefix.ini", argc, argv);
    input.PrintLogo();
    input.PrintParameters();

    // Allocate the grid on device
    Grid grid(input);

    // Allocate the grid on host as a mirror
    GridHost gridHost(grid);

    // Actually make the grid
    gridHost.MakeGrid(input);
    gridHost.SyncToDevice();

    // Make a data array
    DataBlock data;
    data.InitFromGrid(grid, input);

    idfx::cout << "Init Hydrodynamics." << std::endl;
    Hydro hydro(input, grid);

    idfx::cout << "Init Time Integrator." << std::endl;
    TimeIntegrator Tint(input, hydro);

    idfx::cout << "Init Setup." << std::endl;
    Setup mysetup(input,grid,data,hydro);

    idfx::cout << "init Output Routines." << std::endl;
    OutputVTK outVTK(input, data, Tint.getT()); 
    OutputDump outDMP(input, data, Tint.getT());
    
    // Apply initial conditions
    

      // Are we restarting?
    if(input.CheckEntry("CommandLine","restart") > 0) {
      idfx::cout << "Restarting from dump file"  << std::endl;
      outDMP.Read(grid,data,Tint,outVTK,input.GetInt("CommandLine","restart",0));
      hydro.SetBoundary(data,Tint.getT());
      outVTK.Write(data,Tint.getT());
    }
    else {
      idfx::cout << "Creating initial conditions." << std::endl;
      mysetup.InitFlow(data);
      hydro.SetBoundary(data,Tint.getT());
      outDMP.Write(grid, data, Tint, outVTK);
      outVTK.Write(data,Tint.getT());
    }

    idfx::cout << "Cycling Time Integrator..." << std::endl;

    Kokkos::Timer timer;

    real tstop=input.GetReal("TimeIntegrator","tstop",0);

    while(Tint.getT() < tstop) {
      if(tstop-Tint.getT() < Tint.getDt()) Tint.setDt(tstop-Tint.getT());
      Tint.Cycle(data);
      outDMP.Write(grid, data, Tint, outVTK);
      outVTK.Write(data, Tint.getT());
    }
    double tintegration = (timer.seconds()/(grid.np_int[IDIR]*grid.np_int[JDIR]*grid.np_int[KDIR]*Tint.getNcycles()));
    idfx::cout << "Reached t=" << Tint.getT() << std::endl;
    idfx::cout << "Completed in " << timer.seconds() << "seconds and " << Tint.getNcycles() << " cycles. Perfs are " << 1/tintegration << " cell updates/second." << std::endl;
    
    idfx::cout << "Job's done" << std::endl;
  }
  Kokkos::finalize();
#ifdef WITH_MPI    
  MPI_Finalize();
#endif
  return 0;
}
