// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <cstdio>
#include <iomanip>
#include "idefix.hpp"
#include "timeIntegrator.hpp"
#include "input.hpp"


TimeIntegrator::TimeIntegrator(Input & input, DataBlock & data) {
  idfx::pushRegion("TimeIntegrator::TimeIntegrator(Input...)");

  this->timer.reset();
  this->lastLog=timer.seconds();
  this->lastMpiLog=idfx::mpiCallsTimer + idfx::mpiCallsTimer;

  nstages=input.GetInt("TimeIntegrator","nstages",0);

  if(input.CheckEntry("TimeIntegrator","fixed_dt")>0) {
    this->haveFixedDt = true;
    this->fixedDt = input.GetReal("TimeIntegrator","fixed_dt",0);
    idfx::cout << "TimeIntegrator: Using fixed dt time stepping. Ignoring CFL and first_dt."
               << std::endl;
    data.dt=fixedDt;
  }

  if(!haveFixedDt) {
    cfl=input.GetReal("TimeIntegrator","CFL",0);

    if(input.CheckEntry("TimeIntegrator","CFL_max_var")>0) {
      cflMaxVar=input.GetReal("TimeIntegrator","CFL_max_var",0);
    } else {
      cflMaxVar=1.1;
      idfx::cout << "TimeIntegrator: No CFL_max_var defined. Using 1.1 by default." << std::endl;
    }

    if(input.CheckEntry("TimeIntegrator","first_dt")>0) {
      data.dt=input.GetReal("TimeIntegrator","first_dt",0);
    } else {
      data.dt=1e-10;
      idfx::cout << "TimeIntegrator: No first_dt provided. Using dt=1e-10."
                << "This is probably not optimal" << std::endl;
    }
  }

  data.t=0.0;
  ncycles=0;


  if(input.CheckEntry("Output","log")>0) {
    this->cyclePeriod=input.GetInt("Output","log",0);
  }

  if(input.CheckEntry("TimeIntegrator","max_runtime")>0) {
    this->maxRuntime = 3600*input.GetReal("TimeIntegrator","max_runtime",0);
  }

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
  idfx::cout << std::setw(col_width) << data.t;
  idfx::cout << " | " << std::setw(col_width) << ncycles;
  idfx::cout << std::scientific;
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
  idfx::cout << " | " << std::setw(col_width) << divB;
  if(divB>1e-6) {
    IDEFIX_ERROR("TimeIntegrator::Cycle divB>1e-6, check your calculation");
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
  IdefixArray4D<real> Uc = data.hydro.Uc;
  IdefixArray4D<real> Vs = data.hydro.Vs;
  IdefixArray4D<real> Uc0 = data.hydro.Uc0;
  IdefixArray4D<real> Vs0 = data.hydro.Vs0;
  IdefixArray3D<real> InvDt = data.hydro.InvDt;

  real newdt;


  idfx::pushRegion("TimeIntegrator::Cycle");

  //if(timer.seconds()-lastLog >= 1.0) {
  if(ncycles%cyclePeriod==0) ShowLog(data);

  if(haveRKL && (ncycles%2)==1) {    // Runge-Kutta-Legendre cycle
    rkl.Cycle();
  }

  // Apply Boundary conditions
  data.SetBoundaries();

  // Remove Fargo velocity so that the integrator works on the residual
  if(data.haveFargo) data.fargo.SubstractVelocity(data.t);

  // Convert current state into conservative variable and save it
  data.hydro.ConvertPrimToCons();

  // Store initial stage for multi-stage time integrators
  if(nstages>1) {
    Kokkos::deep_copy(Uc0,Uc);
#if MHD == YES
    Kokkos::deep_copy(Vs0,Vs);
#endif
  }

  // save t at the begining of the cycle
  const real t0 = data.t;

  // Reinit datablock for a new stage
  data.ResetStage();

#ifdef WITH_MPI
  MPI_Request dtReduce;
#endif

  for(int stage=0; stage < nstages ; stage++) {
    // If gravity is needed, update it
    if(data.haveGravity) data.gravity.ComputeGravity();

    // Update Uc & Vs
    data.EvolveStage();

    // evolve dt accordingly
    data.t += data.dt;

    // Look for Nans every now and then (this actually cost a lot of time on GPUs
    // because streams are divergent)
    if(ncycles%100==0) if(data.CheckNan()>0) IDEFIX_ERROR("Nan found after integration cycle");

    // Compute next time_step during first stage
    if(stage==0) {
      if(!haveFixedDt) {
        idefix_reduce("Timestep_reduction",
          data.beg[KDIR], data.end[KDIR],
          data.beg[JDIR], data.end[JDIR],
          data.beg[IDIR], data.end[IDIR],
          KOKKOS_LAMBDA (int k, int j, int i, real &dtmin) {
                  dtmin=FMIN(ONE_F/InvDt(k,j,i),dtmin);
              },
          Kokkos::Min<real>(newdt));

        Kokkos::fence();

        newdt=newdt*cfl;

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

      idefix_for("Cycle-update",0,NVAR,0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
        KOKKOS_LAMBDA (int n, int k, int j, int i) {
          Uc(n,k,j,i) = wcs*Uc(n,k,j,i) + w0s*Uc0(n,k,j,i);
      });

#if MHD==YES
      idefix_for("Cycle-update",0,DIMENSIONS,0,data.np_tot[KDIR]+KOFFSET,
                                             0,data.np_tot[JDIR]+JOFFSET,
                                             0,data.np_tot[IDIR]+IOFFSET,
        KOKKOS_LAMBDA (int n, int k, int j, int i) {
          Vs(n,k,j,i) = wcs*Vs(n,k,j,i) + w0s*Vs0(n,k,j,i);
      });
#endif
      // update t
      data.t = wcs*data.t + w0s*t0;

      // Tentatively high order fargo
      //if(data.haveFargo) data.fargo.ShiftSolution(data.t,wcs*data.dt);
    } //else {
      //if(data.haveFargo) data.fargo.ShiftSolution(data.t,data.dt);
    //}

    // Shift solution according to fargo if this is our last stage

    if(data.haveFargo && stage==nstages-1) {
      data.fargo.ShiftSolution(t0,data.dt);
    }


    // Back to using Vc
    data.hydro.ConvertConsToPrim();

    // Check if this is our last stage
    if(stage<nstages-1) {
      // No: Apply boundary conditions & Recompute conservative variables
      // Add back fargo velocity so that boundary conditions are applied on the total V
      if(data.haveFargo) data.fargo.AddVelocity(data.t);
      data.SetBoundaries();
      // And substract it back for next stage
      if(data.haveFargo) data.fargo.SubstractVelocity(data.t);
      data.hydro.ConvertPrimToCons();
    }
  }


  // Add back Fargo velocity so that updated Vc stores the total Velocity
  if(data.haveFargo) data.fargo.AddVelocity(data.t);

  // Wait for hydro/newDt MPI reduction
#ifdef WITH_MPI
  if(idfx::psize>1) {
    MPI_SAFE_CALL(MPI_Wait(&dtReduce, MPI_STATUS_IGNORE));
  }
#endif

  if(haveRKL && (ncycles%2)==0) {    // Runge-Kutta-Legendre cycle
    rkl.Cycle();
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
      IDEFIX_ERROR(msg);
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
  if(this->maxRuntime < 0) return(false);

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
