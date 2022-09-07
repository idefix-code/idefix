// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <cstdio>
#include <iomanip>
#include <string>
#include "idefix.hpp"
#include "timeIntegrator.hpp"
#include "input.hpp"
#include "stateContainer.hpp"


TimeIntegrator::TimeIntegrator(Input & input, DataBlock & data) {
  idfx::pushRegion("TimeIntegrator::TimeIntegrator(Input...)");

  this->timer.reset();
  this->lastLog=timer.seconds();
  this->lastMpiLog=idfx::mpiCallsTimer + idfx::mpiCallsTimer;

  nstages=input.Get<int>("TimeIntegrator","nstages",0);

  if(input.CheckEntry("TimeIntegrator","fixed_dt")>0) {
    this->haveFixedDt = true;
    this->fixedDt = input.Get<real>("TimeIntegrator","fixed_dt",0);
    data.dt=fixedDt;
  }

  if(!haveFixedDt) {
    cfl=input.Get<real>("TimeIntegrator","CFL",0);
    cflMaxVar = input.GetOrSet<real>("TimeIntegrator","CFL_max_var",0, 1.1);
    data.dt = input.GetOrSet<real>("TimeIntegrator","first_dt",0, 1.0e-10);
  }

  this->cyclePeriod = input.GetOrSet<int>("Output","log",0, 100);
  this->maxRuntime = 3600*input.GetOrSet<double>("TimeIntegrator","max_runtime",0.0,-1.0);

  #if defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_CUDA)
    // When Cuda/HIP is enabled, increase the periodicity of nan checks, as this takes a lot of time
    // on GPUs
    checkNanPeriodicity = 100;
  #endif

  // override default check nan periodicity if user decides to
  checkNanPeriodicity = input.GetOrSet<int>("TimeIntegrator","check_nan", 0,
                                            checkNanPeriodicity);

  data.t=0.0;
  ncycles=0;

  if(nstages==2) {
    wc[0] = 0.5;
    w0[0] = 0.5;
  }
  if(nstages==3) {
    wc[0] = 0.25;
    w0[0] = 0.75;
    wc[1] = 2.0/3.0;
    w0[1] = 1.0/3.0;
  }

  // Init the RKL scheme if it's needed
  if(data.hydro.haveRKLParabolicTerms) {
    rkl.Init(input,data);
    haveRKL = true;
  }

  // If multi-stage, create a new state in the datablock called "begin"
  if(nstages>1) {
    data.states["begin"] = StateContainer();
    data.states["begin"].AllocateAs(data.states["current"]);
  }

  idfx::popRegion();
}


void TimeIntegrator::ShowLog(DataBlock &data) {
  if(isSilent) return;
  double rawperf = (timer.seconds()-lastLog)/data.mygrid->np_int[IDIR]/data.mygrid->np_int[JDIR]
                      /data.mygrid->np_int[KDIR]/cyclePeriod;
#ifdef WITH_MPI
  // measure time spent in expensive MPI calls
  double mpiCycleTime = idfx::mpiCallsTimer - lastMpiLog;
  // reduce to an normalized overhead in %
  double mpiOverhead = 100.0 * mpiCycleTime / (timer.seconds() - lastLog);
  lastMpiLog = idfx::mpiCallsTimer;
#endif
  lastLog = timer.seconds();


  int col_width{16};
  if (ncycles == 0) {
    idfx::cout << "TimeIntegrator: ";
    idfx::cout << std::setw(col_width) << "time";
    idfx::cout << " | " << std::setw(col_width) << "cycle";
    idfx::cout << " | " << std::setw(col_width) << "time step";
    idfx::cout << " | " << std::setw(col_width) << "cell (updates/s)";
#ifdef WITH_MPI
    idfx::cout << " | " << std::setw(col_width) << "MPI overhead (%)";
#endif
#if MHD == YES
    idfx::cout << " | " << std::setw(col_width) << "div B";
#endif
    if(haveRKL) {
      idfx::cout << " | " << std::setw(col_width) << "RKL stages";
    }
    idfx::cout << std::endl;
  }

  idfx::cout << "TimeIntegrator: ";
  idfx::cout << std::scientific;
  idfx::cout << std::setw(col_width) << data.t;
  idfx::cout << " | " << std::setw(col_width) << ncycles;
  idfx::cout << " | " << std::setw(col_width) << data.dt;
  if(ncycles>=cyclePeriod) {
    idfx::cout << " | " << std::setw(col_width) << 1 / rawperf;
#ifdef WITH_MPI
  idfx::cout << std::fixed;
    idfx::cout << " | " << std::setw(col_width) << mpiOverhead;
#endif
  } else {
    idfx::cout << " | " << std::setw(col_width) << "N/A";
#if WITH_MPI
    idfx::cout << " | " << std::setw(col_width) << "N/A";
#endif
  }


#if MHD == YES
  // Check divB
  real divB =  data.hydro.CheckDivB();
  idfx::cout << std::scientific;
  idfx::cout << " | " << std::setw(col_width) << divB;
  #ifndef SINGLE_PRECISION
    const real maxdivB = 1e-6;
  #else
    const real maxdivB = 1e-2;
  #endif
  if(divB>maxdivB) {
    std::stringstream msg;
    msg << std::endl << "divB too large, check your calculation";
    throw std::runtime_error(msg.str());
  }
#endif
  if(haveRKL) {
    idfx::cout << " | " << std::setw(col_width) << rkl.stage;
  }
  idfx::cout << std::endl;
}


// Compute one full cycle of the time Integrator
void TimeIntegrator::Cycle(DataBlock &data) {
  // Do one cycle
  IdefixArray3D<real> InvDt = data.hydro.InvDt;
  real newdt;

  idfx::pushRegion("TimeIntegrator::Cycle");

  if(ncycles%cyclePeriod==0) ShowLog(data);

  if(haveRKL && (ncycles%2)==1) {    // Runge-Kutta-Legendre cycle
    rkl.Cycle();
  }

  // save t at the begining of the cycle
  const real t0 = data.t;

  // Reinit datablock for a new stage
  data.ResetStage();

#ifdef WITH_MPI
  MPI_Request dtReduce;
#endif

  /////////////////////////////////////////////////
  // BEGIN STAGES LOOP                           //
  /////////////////////////////////////////////////
  for(int stage=0; stage < nstages ; stage++) {
    // Apply Boundary conditions
    data.SetBoundaries();

    // Remove Fargo velocity so that the integrator works on the residual
    if(data.haveFargo) data.fargo.SubstractVelocity(data.t);

    // Convert current state into conservative variable and save it
    data.hydro.ConvertPrimToCons();

    // Store (deep copy) initial stage for multi-stage time integrators
    if(nstages>1 && stage==0) {
      data.states["begin"].CopyFrom(data.states["current"]);
    }
    // If gravity is needed, update it
    if(data.haveGravity) data.gravity.ComputeGravity();

    // Update Uc & Vs
    data.EvolveStage();

    // evolve dt accordingly
    data.t += data.dt;

    // Look for Nans every now and then (this actually cost a lot of time on GPUs
    // because streams are divergent)
    if(ncycles%checkNanPeriodicity==0) {
      if(data.CheckNan()>0) {
        throw std::runtime_error(std::string("Nan found after integration cycle"));
      }
    }

    // Compute next time_step during first stage
    if(stage==0) {
      if(!haveFixedDt) {
        newdt = cfl*data.ComputeTimestep();
        #ifdef WITH_MPI
          if(idfx::psize>1) {
            MPI_SAFE_CALL(MPI_Iallreduce(MPI_IN_PLACE, &newdt, 1, realMPI, MPI_MIN, MPI_COMM_WORLD,
                                        &dtReduce));
          }
        #endif
      }
    }

    // Is this not the first stage?
    if(stage>0) {
      // do the partial evolution required by the multi-step
      real wcs=wc[stage-1];
      real w0s=w0[stage-1];
      data.states["current"].AddAndStore(wcs, w0s, data.states["begin"]);

      // update t
      data.t = wcs*data.t + w0s*t0;
    }
    // Shift solution according to fargo if this is our last stage
    if(data.haveFargo && stage==nstages-1) {
      data.fargo.ShiftSolution(t0,data.dt);
    }

    // Coarsen conservative variables once they have been evolved
    if(data.haveGridCoarsening) {
      data.Coarsen();
    }

    // Back to using Vc
    data.hydro.ConvertConsToPrim();

    // Add back fargo velocity so that boundary conditions are applied on the total V
    if(data.haveFargo) data.fargo.AddVelocity(data.t);
  }
  /////////////////////////////////////////////////
  // END STAGES LOOP                             //
  /////////////////////////////////////////////////

  // Wait for dt MPI reduction
#ifdef WITH_MPI
  if(!haveFixedDt && idfx::psize>1) {
    MPI_SAFE_CALL(MPI_Wait(&dtReduce, MPI_STATUS_IGNORE));
  }
#endif

  if(haveRKL && (ncycles%2)==0) {    // Runge-Kutta-Legendre cycle
    rkl.Cycle();
  }

  // Coarsen the grid
  if(data.haveGridCoarsening) {
    data.Coarsen();
  }
  // Update current time (should have already been done, but this gets rid of roundoff errors)
  data.t=t0+data.dt;

  if(haveRKL) {
    // update next time step
    real tt = newdt/rkl.dt;
    newdt *= std::fmin(ONE_F, rkl.rmax_par/(tt));
  }

  // Next time step
  if(!haveFixedDt) {
    if(newdt>cflMaxVar*data.dt) {
      data.dt=cflMaxVar*data.dt;
    } else {
      if(ncycles==0 && newdt < 0.5*data.dt) {
        std::stringstream msg;
        msg << "Your guessed first_dt is too large. My next dt=" << newdt << std::endl;
        msg << "Try to reduce first_dt in the ini file.";
        IDEFIX_ERROR(msg);
      }
      data.dt=newdt;
    }
    if(data.dt < 1e-15) {
      std::stringstream msg;
      msg << "dt = " << data.dt << " is too small.";
      throw std::runtime_error(msg.str());
    }
  } else {
    data.dt = fixedDt;
  }


  ncycles++;

  idfx::popRegion();
}



int64_t TimeIntegrator::GetNCycles() {
  return(ncycles);
}

// Check whether our maximumruntime has been reached. Reduce the results on all of the cores
// to make sure they stop simultaneously even if running time are not perfectly in sync
bool TimeIntegrator::CheckForMaxRuntime() {
  idfx::pushRegion("TimeIntegrator::CheckForMaxRuntime");
  // if maxRuntime is negative, this function is disabled (default)
  if(this->maxRuntime < 0) {
    idfx::popRegion();
    return(false);
  }

  double runtime = timer.seconds();
  bool runtimeReached{false};
#ifdef WITH_MPI
  int runtimeValue = 0;
  if(runtime >= this->maxRuntime) runtimeValue = 1;
  MPI_Bcast(&runtimeValue, 1, MPI_INT, 0, MPI_COMM_WORLD);
  runtimeReached = runtimeValue > 0;
#else
  runtimeReached = runtime >= this->maxRuntime;
#endif
  if(runtimeReached) {
    idfx::cout << "TimeIntegrator:CheckForMaxRuntime: Maximum runtime reached."
               << std::endl;
  }
  idfx::popRegion();
  return(runtimeReached);
}

void TimeIntegrator::ShowConfig() {
  if(nstages==1) {
    idfx::cout << "TimeIntegrator: using 1st Order (EULER) integrator." << std::endl;
  } else if(nstages==2) {
    idfx::cout << "TimeIntegrator: using 2nd Order (RK2) integrator." << std::endl;
  } else if(nstages==3) {
    idfx::cout << "TimeIntegrator: using 3rd Order (RK3) integrator." << std::endl;
  } else {
    IDEFIX_ERROR("Unknown time integrator");
  }
  if(haveFixedDt) {
    idfx::cout << "TimeIntegrator: Using fixed dt=" << fixedDt << ". Ignoring CFL and first_dt."
              << std::endl;
  } else {
    idfx::cout << "TimeIntegrator: Using adaptive dt with CFL=" << cfl << " ." << std::endl;
  }
  if(maxRuntime>0) {
    idfx::cout << "TimeIntegrator: will stop after " << maxRuntime/3600 << " hours." << std::endl;
  }

  if(haveRKL) {
    rkl.ShowConfig();
  }
}
