// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

/*
//@HEADER
// ************************************************************************
//
//                        IDEFIX v 1.0-alpha
//
// ************************************************************************
//@HEADER
*/


#include <sys/time.h>

#include <stdlib.h>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <Kokkos_Core.hpp>

#include "idefix.hpp"
#include "profiler.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "gridHost.hpp"
#include "hydro.hpp"
#include "dataBlock.hpp"
#include "timeIntegrator.hpp"
#include "setup.hpp"
#include "output.hpp"
#include "tuner.hpp"




int main( int argc, char* argv[] ) {
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

    Input input = Input(argc, argv);
    input.PrintLogo();
    input.PrintParameters();

    if(initKokkosBeforeMPI) {
      idfx::cout << "Main: detected your configuration needed Kokkos to be initialised before MPI. "
                 << std::endl;
    }

    idfx::cout << "Main: Init Grid." << std::endl;
    // Allocate the grid on device
    Grid grid(input);
    // Allocate the grid image on host
    GridHost gridHost(grid);

    // Actually make the grid on host and sync it on the device
    gridHost.MakeGrid(input);
    gridHost.SyncToDevice();

    // Make a datablock
    idfx::cout << "Main: Init DataBlock." << std::endl;
    DataBlock data;
    data.InitFromGrid(grid, input);

    idfx::cout << "Main: Init Time Integrator." << std::endl;
    TimeIntegrator Tint(input,data);

    idfx::cout << "Main: Init Output Routines." << std::endl;
    Output output(input, data);

    idfx::cout << "Main: Init Setup." << std::endl;
    Setup mysetup(input, grid, data, output);

    // if the user asked for auto-tune, then tune loops now
    if(input.tuningRequested) {
      Tuner::tuneLoops(data, mysetup, input);
    }
    // Apply initial conditions

    // Are we restarting?
    if(input.restartRequested) {
      idfx::cout << "Main: Restarting from dump file."  << std::endl;
      output.RestartFromDump(data,input.restartFileNumber);
      data.SetBoundaries();
      output.CheckForWrites(data);
    } else {
      idfx::cout << "Main: Creating initial conditions." << std::endl;
      idfx::pushRegion("Setup::Initflow");
      mysetup.InitFlow(data);
      idfx::popRegion();
      data.SetBoundaries();
      output.CheckForWrites(data);
      if(data.CheckNan()) {
        IDEFIX_ERROR("Nans were found in your initial conditions.");
      }
    }

    idfx::cout << "Main: Cycling Time Integrator..." << std::endl;

    Kokkos::Timer timer;
    output.ResetTimer();

    real tstop = input.GetReal("TimeIntegrator","tstop",0);

    while(data.t < tstop) {
      if(tstop-data.t < data.dt) data.dt = tstop-data.t;
      Tint.Cycle(data);
      output.CheckForWrites(data);
      if(input.CheckForAbort() || Tint.CheckForMaxRuntime() ) {
        idfx::cout << "Main: Saving current state and aborting calculation." << std::endl;
        output.ForceWrite(data);
        break;
      }
      if(input.maxCycles>=0) {
        if(Tint.GetNCycles() >= input.maxCycles) {
          idfx::cout << "Main: Reached maximum number of integration cycles." << std::endl;
          break;
        }
      }
    }

    int n_days{0}, n_hours{0}, n_minutes{0}, n_seconds{0};
    div_t divres;
    divres = div(timer.seconds(), 86400);
    n_days = divres.quot;
    divres = div(divres.rem, 3600);
    n_hours = divres.quot;
    divres = div(divres.rem, 60);
    n_minutes = divres.quot;
    n_seconds = divres.rem;

    double tintegration = timer.seconds() / grid.np_int[IDIR] / grid.np_int[JDIR]
                            / grid.np_int[KDIR] / Tint.GetNCycles();

    idfx::cout << "Main: Reached t=" << data.t << std::endl;
    idfx::cout << "Main: Completed in ";
    if (n_days > 0) {
      idfx::cout << n_days << " day";
      if (n_days != 1) {
        idfx::cout << "s";
      }
      idfx::cout << " ";
    }
    if (n_hours > 0) {
      idfx::cout << n_hours << " hour";
      if (n_hours != 1) {
        idfx::cout << "s";
      }
      idfx::cout << " ";
    }
    if (n_minutes > 0) {
      idfx::cout << n_minutes << " minute";
      if (n_minutes != 1) {
        idfx::cout << "s";
      }
      idfx::cout << " ";
    }
    idfx::cout << n_seconds << " second";
    if (n_seconds != 1) {
      idfx::cout << "s";
    }
    idfx::cout << " ";
    idfx::cout << "and " << Tint.GetNCycles() << " cycle";
    if (Tint.GetNCycles() != 1) {
      idfx::cout << "s";
    }
    idfx::cout << std::endl;
    idfx::cout << "Main: ";
    idfx::cout << "Perfs are " << 1/tintegration << " cell updates/second" << std::endl;
    #ifdef WITH_MPI
      idfx::cout << "MPI overhead represents "
                 << static_cast<int>(100.0*idfx::mpiCallsTimer/timer.seconds())
                 << "% of total run time." << std::endl;
    #endif

    idfx::cout << "Outputs represent "
               << static_cast<int>(100.0*output.GetTimer()/timer.seconds())
              << "% of total run time." << std::endl;
    // Show profiler output
    idfx::prof.Show();
  }
  idfx::cout << "Main: Job's done" << std::endl;
  Kokkos::finalize();

#ifdef WITH_MPI
  MPI_Finalize();
#endif

  return 0;
}
