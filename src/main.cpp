// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

/*
//@HEADER
// ************************************************************************
//
//                        IDEFIX v 0.0-alpha
//
// ************************************************************************
//@HEADER
*/


#include <sys/time.h>

#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <csignal>

#include <Kokkos_Core.hpp>

#include "idefix.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "gridHost.hpp"
#include "hydro.hpp"
#include "dataBlock.hpp"
#include "timeIntegrator.hpp"
#include "setup.hpp"
#include "output.hpp"


bool abortRequested;

void signalHandler(int signum) {
  idfx::cout << std::endl << "Main:: Caught interrupt " << signum << std::endl;
  abortRequested=true;
}


int main( int argc, char* argv[] ) {
  bool haveSlurm = false;

  // When running on GPUS allocated from a SLURM scheduler,
  // Kokkos needs to be initialised *before* the MPI layer
#ifdef KOKKOS_ENABLE_CUDA
  if(std::getenv("SLURM_JOB_ID") != NULL) {
    haveSlurm = true;
  }
#endif

  if(haveSlurm)  Kokkos::initialize( argc, argv );

#ifdef WITH_MPI
  MPI_Init(&argc,&argv);
#endif

  if(!haveSlurm) Kokkos::initialize( argc, argv );


  {
    idfx::initialize();

    //signal(SIGINT, signalHandler);
    //signal(SIGTERM, signalHandler);
    signal(SIGUSR2, signalHandler);
    abortRequested=false;

    Input input = Input(argc, argv);
    input.PrintLogo();
    input.PrintParameters();

    if(haveSlurm)
      idfx::cout << "Main:: detected you were runnning on a SLURM scheduler with CUDA. "
                    "Kokkos has been initialised accordingly." << std::endl;

    idfx::cout << "Main::Init Grid." << std::endl;
    // Allocate the grid on device
    Grid grid(input);
    // Allocate the grid image on host
    GridHost gridHost(grid);

    // Actually make the grid on host and sync it on the device
    gridHost.MakeGrid(input);
    gridHost.SyncToDevice();

    // Make a datablock
    idfx::cout << "Main::Init DataBlock." << std::endl;
    DataBlock data;
    data.InitFromGrid(grid, input);

    idfx::cout << "Main::Init Time Integrator." << std::endl;
    TimeIntegrator Tint(input,data);

    idfx::cout << "Main::Init Setup." << std::endl;
    Setup mysetup(input,grid,data);

    idfx::cout << "Main::Onit Output Routines." << std::endl;
    Output output(input, data, mysetup);

    // Apply initial conditions

    // Are we restarting?
    if(input.CheckEntry("CommandLine","restart") > 0) {
      idfx::cout << "Main::Restarting from dump file"  << std::endl;
      output.RestartFromDump(data,input.GetInt("CommandLine","restart",0));
      data.hydro.SetBoundary(data.t);
      output.CheckForWrites(data);
    } else {
      idfx::cout << "Main::Creating initial conditions." << std::endl;
      mysetup.InitFlow(data);
      data.hydro.SetBoundary(data.t);
      output.CheckForWrites(data);
    }

    idfx::cout << "Main::Cycling Time Integrator..." << std::endl;

    Kokkos::Timer timer;

    real tstop = input.GetReal("TimeIntegrator","tstop",0);

    while(data.t < tstop) {
      if(tstop-data.t < data.dt) data.dt = tstop-data.t;
      Tint.Cycle(data);
      output.CheckForWrites(data);
      if(abortRequested) {
        idfx::cout << "Main:: Saving current state and aborting calculation" << std::endl;
        output.ForceWrite(data);
        break;
      }
    }

    double tintegration = timer.seconds() / grid.np_int[IDIR] / grid.np_int[JDIR]
                            / grid.np_int[KDIR] / Tint.getNcycles();

    idfx::cout << "Main::Reached t=" << data.t << std::endl;
    idfx::cout << "Main::Completed in " << timer.seconds() << "seconds and " << Tint.getNcycles()
               << " cycles. Perfs are " << 1/tintegration << " cell updates/second." << std::endl;

    idfx::cout << "Main::Job's done" << std::endl;
  }

  Kokkos::finalize();

#ifdef WITH_MPI
  MPI_Finalize();
#endif

  return 0;
}
