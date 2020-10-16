// ********************************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ********************************************************************************************************

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
#include <csignal>

#include <Kokkos_Core.hpp>

#include "idefix.hpp"

bool abortRequested;

void signalHandler(int signum) {
  idfx::cout << std::endl << "Main:: Caught interrupt " << signum << std::endl;
  abortRequested=true;
}


int main( int argc, char* argv[] )
{

// When running on GPUS allocated from a SLURM scheduler, Kokkos needs to be initialised *before* the MPI layer
#ifdef SLURM_MPI
    Kokkos::initialize( argc, argv );
#endif

#ifdef WITH_MPI
  MPI_Init(&argc,&argv);
#endif

#ifndef SLURM_MPI
  Kokkos::initialize( argc, argv );
#endif

  {

    idfx::initialize();

    //signal(SIGINT, signalHandler);
    //signal(SIGTERM, signalHandler);
    signal(SIGUSR2, signalHandler);
    abortRequested=false;

    Input input = Input(argc, argv);
    input.PrintLogo();
    input.PrintParameters();

    // Allocate the grid on device
    Grid grid(input);

    // Allocate the grid on host as a mirror
    GridHost gridHost(grid);

    // Actually make the grid
    gridHost.MakeGrid(input);
    gridHost.SyncToDevice();

    idfx::cout << "Main::Init Hydrodynamics." << std::endl;
    Hydro hydro(input, grid);

    idfx::cout << "Main::Init Time Integrator." << std::endl;
    TimeIntegrator Tint(input, hydro);

    // Make a datablock
    DataBlock data;
    data.InitFromGrid(grid, hydro, input);

    idfx::cout << "Main::Init Setup." << std::endl;
    Setup mysetup(input,grid,data,hydro);

    idfx::cout << "Main::Onit Output Routines." << std::endl;
    OutputVTK outVTK(input, data, Tint.getT());
    OutputDump outDMP(input, data, Tint.getT());

    // Apply initial conditions


      // Are we restarting?
    if(input.CheckEntry("CommandLine","restart") > 0) {
      idfx::cout << "Main::Restarting from dump file"  << std::endl;
      outDMP.Read(grid,data,Tint,outVTK,input.GetInt("CommandLine","restart",0));
      hydro.SetBoundary(data,Tint.getT());
      outVTK.CheckForWrite(data,Tint.getT());
    }
    else {
      idfx::cout << "Main::Creating initial conditions." << std::endl;
      mysetup.InitFlow(data);
      hydro.SetBoundary(data,Tint.getT());
      outDMP.CheckForWrite(grid, data, Tint, outVTK);
      outVTK.CheckForWrite(data,Tint.getT());
    }

    idfx::cout << "Main::Cycling Time Integrator..." << std::endl;

    Kokkos::Timer timer;

    real tstop=input.GetReal("TimeIntegrator","tstop",0);

    while(Tint.getT() < tstop) {
      if(tstop-Tint.getT() < Tint.getDt()) Tint.setDt(tstop-Tint.getT());
      Tint.Cycle(data);
      outDMP.CheckForWrite(grid, data, Tint, outVTK);
      outVTK.CheckForWrite(data, Tint.getT());
      if(abortRequested) {
        idfx::cout << "Main:: Saving current state and aborting calculation" << std::endl;
        outDMP.Write(grid, data, Tint, outVTK);
        break;
      }
    }
    double tintegration = timer.seconds()/grid.np_int[IDIR]/grid.np_int[JDIR]/grid.np_int[KDIR]/Tint.getNcycles();
    idfx::cout << "Main::Reached t=" << Tint.getT() << std::endl;
    idfx::cout << "Main::Completed in " << timer.seconds() << "seconds and " << Tint.getNcycles() << " cycles. Perfs are " << 1/tintegration << " cell updates/second." << std::endl;

    idfx::cout << "Main::Job's done" << std::endl;
  }
  Kokkos::finalize();
#ifdef WITH_MPI
  MPI_Finalize();
#endif
  return 0;
}
