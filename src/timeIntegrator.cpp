/// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <cstdio>
#include "idefix.hpp"
#include "timeIntegrator.hpp"

TimeIntegrator::TimeIntegrator(Input & input, Hydro &physics) {
  idfx::pushRegion("TimeIntegrator::TimeIntegrator(Input...)");

  this->hydro=&physics;
  this->timer.reset();
  this->lastLog=timer.seconds();
  this->lastMpiLog=idfx::mpiTimer;

  nstages=input.GetInt("TimeIntegrator","nstages",0);

  dt=input.GetReal("TimeIntegrator","first_dt",0);

  t=0.0;
  ncycles=0;
  cfl=input.GetReal("TimeIntegrator","CFL",0);

  if(input.CheckEntry("Output","log")>0) {
    this->cyclePeriod=input.GetInt("Output","log",0);
  } else {
    this->cyclePeriod = 100; // Default log every 100 loops
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

Hydro& TimeIntegrator::GetHydro() {
  return (*this->hydro);
}

// Compute one Stage of the time Integrator
void TimeIntegrator::Stage(DataBlock &data) {
  idfx::pushRegion("TimeIntegrator::Stage");


  // Compute current when needed
  if(hydro->needCurrent) hydro->CalcCurrent(data);

  // Loop on all of the directions
  for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
    // Step one: extrapolate the variables to the sides, result is stored in the physics object
    hydro->ExtrapolatePrimVar(data, dir);

    // Step 2: compute the intercell flux with our Riemann solver, store the resulting InvDt
    hydro->CalcRiemannFlux(data, dir, t);


    // Step 2.5: compute intercell parabolic flux when needed
    if(hydro->haveParabolicTerms) hydro->CalcParabolicFlux(data, dir, t);

    // Step 3: compute the resulting evolution of the conserved variables, stored in Uc
    hydro->CalcRightHandSide(data, dir, t, dt);
  }

  // Step 4: add source terms to the conserved variables (curvature, rotation, etc)
  if(hydro->haveSourceTerms) hydro->AddSourceTerms(data, t, dt);

#if MHD == YES && DIMENSIONS >= 2
  // Compute the field evolution according to CT
  hydro->CalcCornerEMF(data, t);
  if(hydro->haveResistivity || hydro->haveAmbipolar) hydro->CalcNonidealEMF(data,t);
  hydro->EvolveMagField(data, t, dt);
#endif

  idfx::popRegion();
}

void TimeIntegrator::ReinitInvDt(DataBlock & data) {
  idfx::pushRegion("TimeIntegrator::ReinitInvDt");

  IdefixArray3D<real> InvDt=data.InvDt;

  idefix_for("InitInvDt",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
          InvDt(k,j,i) = ZERO_F;
      });

  idfx::popRegion();
}

// Compute one full cycle of the time Integrator
void TimeIntegrator::Cycle(DataBlock & data) {
  // Do one cycle
  IdefixArray4D<real> Uc = data.Uc;
  IdefixArray4D<real> Vs = data.Vs;
  IdefixArray4D<real> Uc0 = data.Uc0;
  IdefixArray4D<real> Vs0 = data.Vs0;
  IdefixArray3D<real> InvDt=data.InvDt;

  real newdt;

  idfx::pushRegion("TimeIntegrator::Cycle");

  //if(timer.seconds()-lastLog >= 1.0) {
  if(ncycles%cyclePeriod==0) {
    double rawperf = (timer.seconds()-lastLog)/(data.mygrid->np_int[IDIR]*data.mygrid->np_int[JDIR]
                      *data.mygrid->np_int[KDIR]*cyclePeriod);
    #ifdef WITH_MPI
      // measure the time spent in the MPI calls
      double mpiOverhead = (idfx::mpiTimer-lastMpiLog)/
                            (timer.seconds()-lastLog-idfx::mpiTimer+lastMpiLog)*100;
      lastMpiLog = idfx::mpiTimer;
    #endif
    lastLog = timer.seconds();

    idfx::cout << "TimeIntegrator: t=" << t << " Cycle " << ncycles << " dt=" << dt << std::endl;
    if(ncycles>=cyclePeriod) {
      idfx::cout << "\t " << 1/rawperf << " cell updates/second";
      #ifdef WITH_MPI
        idfx::cout << " ; " << mpiOverhead << "% MPI overhead";
      #endif
      idfx::cout << std::endl;
    }
    #if MHD == YES
    // Check divB
    idfx::cout << "\t maxdivB=" << hydro->CheckDivB(data) << std::endl;
    #endif
  }

    // Apply Boundary conditions
  hydro->SetBoundary(data,t);

  // Convert current state into conservative variable and save it
  hydro->ConvertPrimToCons(data);


  // Store initial stage for multi-stage time integrators
  if(nstages>1) {
    Kokkos::deep_copy(Uc0,Uc);
  #if MHD == YES
    Kokkos::deep_copy(Vs0,Vs);
  #endif
  }

  // Reinit timestep
  ReinitInvDt(data);

  for(int stage=0; stage < nstages ; stage++) {
    // Update Uc & Vs
    Stage(data);
    #if MHD == YES && DIMENSIONS >= 2
    hydro->ReconstructVcField(data, data.Uc);
    #endif

    // Look for Nans every now and then (this actually cost a lot of time on GPUs
    // because streams are divergent)
    if(ncycles%100==0) if(data.CheckNan()>0) IDEFIX_ERROR("Nan found after integration cycle");

    // Compute next time_step during first stage
    if(stage==0) {
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
    }
    hydro->ConvertConsToPrim(data);
    
    // Check if this is our last stage
    if(stage<nstages-1) {
      // No: Apply boundary conditions & Recompute conservative variables
        hydro->SetBoundary(data,t);
        hydro->ConvertPrimToCons(data);
    }

  }

  // Update current time
  t=t+dt;


  // Next time step
  if(newdt>1.1*dt) {
    dt=1.1*dt;
  } else {
    dt=newdt;
  }

  ncycles++;

  idfx::popRegion();
}

real TimeIntegrator::getDt() {
  return(dt);
}

real TimeIntegrator::getT() {
  return (t);
}

void TimeIntegrator::setDt(real dtin) {
  dt=dtin;
}

int64_t TimeIntegrator::getNcycles() {
  return(ncycles);
}
