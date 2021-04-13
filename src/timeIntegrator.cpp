// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
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
  this->lastMpiLog=idfx::mpiTimer;

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

  idfx::popRegion();
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
  if(ncycles%cyclePeriod==0) {
    double rawperf = (timer.seconds()-lastLog)/data.mygrid->np_int[IDIR]/data.mygrid->np_int[JDIR]
                      /data.mygrid->np_int[KDIR]/cyclePeriod;
#ifdef WITH_MPI
    // measure the time spent in the MPI calls
    double mpiOverhead = (idfx::mpiTimer-lastMpiLog)
                            / (timer.seconds()-lastLog-idfx::mpiTimer+lastMpiLog)*100;
    lastMpiLog = idfx::mpiTimer;
#endif
    lastLog = timer.seconds();


    int col_width{16};
    if (ncycles == 0) {
      idfx::cout << "TimeIntegrator: ";
      idfx::cout << std::setw(col_width) << "time";
      idfx::cout << " | " << std::setw(col_width) << "cycle";
      idfx::cout << " | " << std::setw(col_width) << "time step";
      idfx::cout << " | " << std::setw(col_width) << "cell updates/s";
#ifdef WITH_MPI
      idfx::cout << " | " << std::setw(col_width) << "MPI overhead (%)";
#endif
#if MHD == YES
      idfx::cout << " | " << std::setw(col_width) << "div B";
#endif
      idfx::cout << std::endl;
    }

    idfx::cout << "TimeIntegrator: ";
    idfx::cout << std::setw(col_width) << data.t;
    idfx::cout << " | " << std::setw(col_width) << ncycles;
    idfx::cout << " | " << std::setw(col_width) << data.dt;
    if(ncycles>=cyclePeriod) {
      idfx::cout << " | " << std::setw(col_width) << 1 / rawperf;
#ifdef WITH_MPI
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
    if(divB>1e-10) {
      IDEFIX_ERROR("TimeIntegrator::Cycle divB>1e-10, check your calculation");
    }
#endif
    idfx::cout << std::endl;
  }

    // Apply Boundary conditions
    // TODO(lesurg): Make a general boundary condition call in datablock
  data.hydro.SetBoundary(data.t);

  // Remove Fargo velocity so that the integrator works on the residual
  if(data.hydro.haveFargo) data.hydro.fargo.SubstractVelocity(data.t);

  // Convert current state into conservative variable and save it
  data.hydro.ConvertPrimToCons();



  // Store initial stage for multi-stage time integrators
  if(nstages>1) {
    Kokkos::deep_copy(Uc0,Uc);
#if MHD == YES
    Kokkos::deep_copy(Vs0,Vs);
#endif
  }

  // Reinit datablock for a new stage
  data.ResetStage();

  for(int stage=0; stage < nstages ; stage++) {
    // Update Uc & Vs
    data.EvolveStage();

    // Look for Nans every now and then (this actually cost a lot of time on GPUs
    // because streams are divergent)
    if(ncycles%100==0) if(data.CheckNan()>0) IDEFIX_ERROR("Nan found after integration cycle");

    // Compute next time_step during first stage
    if(stage==0) {
      if(!haveFixedDt) {
        Kokkos::parallel_reduce("Timestep_reduction",
            Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>
            ({0,0,0},{data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]}),
            KOKKOS_LAMBDA (int k, int j, int i, real &dtmin) {
                dtmin=FMIN(ONE_F/InvDt(k,j,i),dtmin);
            }, Kokkos::Min<real>(newdt) );

        Kokkos::fence();

        newdt=newdt*cfl;

  #ifdef WITH_MPI
        if(idfx::psize>1) {
          MPI_SAFE_CALL(MPI_Allreduce(MPI_IN_PLACE, &newdt, 1, realMPI, MPI_MIN, MPI_COMM_WORLD));
        }
  #endif
      }
    }

    // Is this not the first stage?
    if(stage>0) {
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
      // Tentatively high order fargo
      //if(data.hydro.haveFargo) data.hydro.fargo.ShiftSolution(data.t,wcs*data.dt);
    } //else {
      //if(data.hydro.haveFargo) data.hydro.fargo.ShiftSolution(data.t,data.dt);
    //}

    // Shift solution according to fargo if this is our last stage

    if(data.hydro.haveFargo && stage==nstages-1) {
      data.hydro.fargo.ShiftSolution(data.t,data.dt);
    }


    // Back to using Vc
    data.hydro.ConvertConsToPrim();

    // Check if this is our last stage
    if(stage<nstages-1) {
      // No: Apply boundary conditions & Recompute conservative variables
      // Add back fargo velocity so that boundary conditions are applied on the total V
      if(data.hydro.haveFargo) data.hydro.fargo.AddVelocity(data.t);
      data.hydro.SetBoundary(data.t);
      // And substract it back for next stage
      if(data.hydro.haveFargo) data.hydro.fargo.SubstractVelocity(data.t);
      data.hydro.ConvertPrimToCons();
    }
  }


  // Add back Fargo velocity so that updated Vc stores the total Velocity
  if(data.hydro.haveFargo) data.hydro.fargo.AddVelocity(data.t);


  // Runge-Kutta-Legendre cycle
  data.rkl.Cycle();


  // Update current time
  data.t=data.t+data.dt;

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



int64_t TimeIntegrator::getNcycles() {
  return(ncycles);
}

// Check whether our maximumruntime has been reached. Reduce the results on all of the cores
// to make sure they stop simultaneously even if running time are not perfectly in sync
bool TimeIntegrator::CheckForMaxRuntime() {
  // if maxRuntime is negative, this function is disabled (default)
  if(this->maxRuntime < 0) return(false);

  double runtime = timer.seconds();
  bool runtimeReached{false};
#ifdef WITH_MPI
  int runtimeValue = 0;
  if(runtime >= this->maxRuntime) runtimeValue = 1;
  MPI_Allreduce(MPI_IN_PLACE, &runtimeValue, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  runtimeReached = runtimeValue > 0;
#else
  runtimeReached = runtime >= this->maxRuntime;
#endif
  if(runtimeReached) {
    idfx::cout << "TimeIntegrator:CheckForMaxRuntime: Maximum runtime reached."
               << std::endl;
  }
  return(runtimeReached);
}
