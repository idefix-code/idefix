// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <cmath>

#include "rkl.hpp"
#include "dataBlock.hpp"
#include "hydro.hpp"


#ifndef RKL_ORDER
  #define RKL_ORDER       2
#endif
#define NVAR_MAX         10


void RKLegendre::AddVariable(int var, IdefixArray1D<int>::HostMirror &varListHost ) {
  bool haveit{false};

  // Check whether we have this variable in the list
  for(int i = 0 ; i < nvarRKL ; i++) {
    if(varListHost(i) == var) haveit=true;
  }

  // We don't have it, then add it to the list
  if(!haveit) {
    if(nvarRKL >= NVAR_MAX) IDEFIX_ERROR("RKL: cannot compute that many variables");
    varListHost(nvarRKL) = var;
    nvarRKL++;
  }
}

// Copy just the variables required by the RK scheme
void RKLegendre::Copy(IdefixArray4D<real> &out, IdefixArray4D<real> &in) {
  IdefixArray1D<int> vars = this->varList;

  idefix_for("RKL_Copy",
             0, nvarRKL,
             0, data->np_tot[KDIR],
             0, data->np_tot[JDIR],
             0, data->np_tot[IDIR],
             KOKKOS_LAMBDA(int n, int k, int j, int i) {
               const int var = vars(n);
               out(var,k,j,i) = in(var,k,j,i);
             });
}


void RKLegendre::Init(Input &input, DataBlock &datain) {
  idfx::pushRegion("RKLegendre::Init");

  idfx::cout << "RKLegendre: enabled." << std::endl;
  // Save the datablock to which we are attached from now on
  this->data = &datain;


  if(input.CheckEntry("RKL","cfl")>0) {
    cfl_rkl = input.GetReal("RKL","cfl",0);
  } else {
    idfx::cout << "RKLegendre: no RKL cfl given. Using 0.5 by default." << std::endl;
    cfl_rkl = 0.5;
  }
  rmax_par = 100.0;

  // Make a list of variables
  varList = IdefixArray1D<int>("RKL_VarList",NVAR_MAX);
  IdefixArray1D<int>::HostMirror varListHost = Kokkos::create_mirror_view(varList);

  // Create a list of variables
  // Viscosity
  if(data->hydro.viscosityStatus.isRKL) {
    EXPAND( AddVariable(MX1, varListHost);   ,
            AddVariable(MX2, varListHost);   ,
            AddVariable(MX3, varListHost);   )

    #if HAVE_ENERGY
      AddVariable(ENG, varListHost);
    #endif
  }
  // Ambipolar diffusion
  if(data->hydro.ambipolarStatus.isRKL || data->hydro.resistivityStatus.isRKL) {
    #if COMPONENTS == 3 && DIMENSIONS < 3
      AddVariable(BX3, varListHost);
    #endif
    #if COMPONENTS >= 2 && DIMENSIONS < 2
      AddVariable(BX2, varListHost);
    #endif
    #if HAVE_ENERGY
      AddVariable(ENG, varListHost);
    #endif
    haveVs = true;
  }

  // Copy the list on the device
  Kokkos::deep_copy(varList,varListHost);

  #ifdef WITH_MPI
    if(haveVs) {
      mpi.Init(&datain, data->hydro.Vc, varListHost, nvarRKL, haveVs, data->hydro.Vs);
    } else {
      mpi.Init(&datain, data->hydro.Vc, varListHost, nvarRKL);
    }
  #endif


  // Variable allocation

  dU = IdefixArray4D<real>("RKL_dU", NVAR,
                           data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  dU0 = IdefixArray4D<real>("RKL_dU0", NVAR,
                           data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  Uc1 = IdefixArray4D<real>("RKL_Uc1", NVAR,
                           data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

  if(haveVs) {
    dB = IdefixArray4D<real>("RKL_dB", DIMENSIONS,
                      data->np_tot[KDIR]+KOFFSET,
                      data->np_tot[JDIR]+JOFFSET,
                      data->np_tot[IDIR]+IOFFSET);
    dB0 = IdefixArray4D<real>("RKL_dB0", DIMENSIONS,
                      data->np_tot[KDIR]+KOFFSET,
                      data->np_tot[JDIR]+JOFFSET,
                      data->np_tot[IDIR]+IOFFSET);
    Vs1 = IdefixArray4D<real>("RKL_Vs1", DIMENSIONS,
                      data->np_tot[KDIR]+KOFFSET,
                      data->np_tot[JDIR]+JOFFSET,
                      data->np_tot[IDIR]+IOFFSET);
  }

  idfx::popRegion();
}


void RKLegendre::Cycle() {
  idfx::pushRegion("RKLegendre::Cycle");

  IdefixArray4D<real> dU = this->dU;
  IdefixArray4D<real> dU0 = this->dU0;
  IdefixArray4D<real> Uc = data->hydro.Uc;
  IdefixArray4D<real> Uc0 = data->hydro.Uc0;
  IdefixArray4D<real> Uc1 = this->Uc1;

  IdefixArray4D<real> dB = this->dB;
  IdefixArray4D<real> dB0 = this->dB0;
  IdefixArray4D<real> Vs = data->hydro.Vs;
  IdefixArray4D<real> Vs0 = data->hydro.Vs0;
  IdefixArray4D<real> Vs1 = this->Vs1;

  IdefixArray1D<int> varList = this->varList;
  real time = data->t;

  real dt_hyp = data->dt;

  // Tell the datablock that we're performing the RKL cycle
  data->rklCycle = true;

  // first RKL stage
  stage = 1;

  // Apply Boundary conditions on the full set of variables
  data->hydro.boundary.SetBoundaries(time);

  // Convert current state into conservative variable
  data->hydro.ConvertPrimToCons();

  // Store the result in Uc0
  Copy(Uc0,Uc);
  if(haveVs) Kokkos::deep_copy(Vs0,Vs);

  // evolve RKL stage
  EvolveStage(time);

  ComputeDt();

  Copy(dU0,dU);
  if(haveVs) Kokkos::deep_copy(dB0,dB);

  // Compute number of RKL steps
  real nrkl;
  real scrh =  dt_hyp/dt;
#if RKL_ORDER == 1
  // Solution of quadratic Eq.
  // 2*dt_hyp/dt_exp = s^2 + s
  nrkl = 4.0*scrh / (1.0 + std::sqrt(1.0 + 8.0*scrh));
#elif RKL_ORDER == 2
  // Solution of quadratic Eq.
  // 4*dt_hyp/dt_exp = s^2 + s - 2
  nrkl = 4.0*(1.0 + 2.0*scrh)
          / (1.0 + sqrt(9.0 + 16.0*scrh));
#else
  //#error Invalid RKL_ORDER
#endif
  int rklstages = 1 + floor(nrkl);

  // Compute coefficients
  real w1, mu_tilde_j;
#if RKL_ORDER == 1
  w1 = 2.0/(rklstages*rklstages + rklstages);
  mu_tilde_j = w1;
#elif RKL_ORDER == 2
  real b_j, b_jm1, b_jm2, a_jm1;
  w1 = 4.0/(rklstages*rklstages + rklstages - 2.0);
  mu_tilde_j = w1/3.0;

  b_j = b_jm1 = b_jm2 = 1.0/3.0;
  a_jm1 = 1.0 - b_jm1;
#endif

#if RKL_ORDER == 1
  time = data->t + 0.5*dt_hyp*(stage*stage+stage)*w1;
#elif RKL_ORDER == 2
  time = data->t + 0.25*dt_hyp*(stage*stage+stage-2)*w1;
#endif

  idefix_for("RKL_Cycle_InitUc1",
             0, nvarRKL,
             data->beg[KDIR],data->end[KDIR],
             data->beg[JDIR],data->end[JDIR],
             data->beg[IDIR],data->end[IDIR],
    KOKKOS_LAMBDA (int n, int k, int j, int i) {
      int nv = varList(n);
      Uc1(nv,k,j,i) = Uc(nv,k,j,i);
      Uc(nv,k,j,i) = Uc1(nv,k,j,i) + mu_tilde_j*dt_hyp*dU0(nv,k,j,i);
    }
  );
  if(haveVs) {
    idefix_for("RKL_Cycle_InitVs1",
             0, DIMENSIONS,
             data->beg[KDIR],data->end[KDIR]+KOFFSET,
             data->beg[JDIR],data->end[JDIR]+JOFFSET,
             data->beg[IDIR],data->end[IDIR]+IOFFSET,
    KOKKOS_LAMBDA (int n, int k, int j, int i) {
      Vs1(n,k,j,i) = Vs(n,k,j,i);
      Vs(n,k,j,i) = Vs1(n,k,j,i) + mu_tilde_j*dt_hyp*dB0(n,k,j,i);
    });
    data->hydro.boundary.ReconstructVcField(Uc);
  }

  // Convert current state into primitive variable
  data->hydro.ConvertConsToPrim();

  real mu_j, nu_j, gamma_j;
  // subStages loop
  for(stage=2; stage <= rklstages ; stage++) {
    //idfx::cout << "RKL: looping stages" << std::endl;
    // compute RKL coefficients
#if RKL_ORDER == 1
    mu_j       = (2.0*stage -1.0)/stage;
    mu_tilde_j = w1*mu_j;
    nu_j       = -(stage -1.0)/stage;
#elif RKL_ORDER == 2
    mu_j       = (2.0*stage -1.0)/stage * b_j/b_jm1;
    mu_tilde_j = w1*mu_j;
    gamma_j    = -a_jm1*mu_tilde_j;
    nu_j       = -(stage -1.0)*b_j/(stage*b_jm2);

    b_jm2 = b_jm1;
    b_jm1 = b_j;
    a_jm1 = 1.0 - b_jm1;
    b_j   = 0.5*(stage*stage+3.0*stage)/(stage*stage+3.0*stage+2.0);
#endif

    // Apply Boundary conditions
    this->SetBoundaries(time);

    // evolve RKL stage
    EvolveStage(time);

    // update Uc
    idefix_for("RKL_Cycle_UpdateUc",
             0, nvarRKL,
             data->beg[KDIR],data->end[KDIR],
             data->beg[JDIR],data->end[JDIR],
             data->beg[IDIR],data->end[IDIR],
      KOKKOS_LAMBDA (int n, int k, int j, int i) {
        const int nv = varList(n);
        real Y = mu_j*Uc(nv,k,j,i) + nu_j*Uc1(nv,k,j,i);
        Uc1(nv,k,j,i) = Uc(nv,k,j,i);
#if RKL_ORDER == 1
        Uc(nv,k,j,i) = Y + dt_hyp*mu_tilde_j*dU(nv,k,j,i);
#elif RKL_ORDER == 2
        Uc(nv,k,j,i) = Y + (1.0 - mu_j - nu_j)*Uc0(nv,k,j,i)
                                + dt_hyp*mu_tilde_j*dU(nv,k,j,i)
                                + gamma_j*dt_hyp*dU0(nv,k,j,i);
#endif
        });

    if(haveVs) {
      // update Vs
      idefix_for("RKL_Cycle_UpdateVs",
              0, DIMENSIONS,
              data->beg[KDIR],data->end[KDIR]+KOFFSET,
              data->beg[JDIR],data->end[JDIR]+JOFFSET,
              data->beg[IDIR],data->end[IDIR]+IOFFSET,
        KOKKOS_LAMBDA (int n, int k, int j, int i) {
          real Y = mu_j*Vs(n,k,j,i) + nu_j*Vs1(n,k,j,i);
          Vs1(n,k,j,i) = Vs(n,k,j,i);
  #if RKL_ORDER == 1
          Vs(n,k,j,i) = Y + dt_hyp*mu_tilde_j*dB(n,k,j,i);
  #elif RKL_ORDER == 2
          Vs(n,k,j,i) = Y + (1.0 - mu_j - nu_j)*Vs0(n,k,j,i)
                                  + dt_hyp*mu_tilde_j*dB(n,k,j,i)
                                  + gamma_j*dt_hyp*dB0(n,k,j,i);
  #endif
          });
      data->hydro.boundary.ReconstructVcField(Uc);
    }
    // Convert current state into primitive variable
    data->hydro.ConvertConsToPrim();

    // increment time
#if RKL_ORDER == 1
    time = data->t + 0.5*dt_hyp*(stage*stage+stage)*w1;
#elif RKL_ORDER == 2
    time = data->t + 0.25*dt_hyp*(stage*stage+stage-2)*w1;
#endif
  }

  // Tell the datablock that we're done
  data->rklCycle = false;
  idfx::popRegion();
}


void RKLegendre::ResetFlux() {
  idfx::pushRegion("RKLegendre::ResetFlux");
  IdefixArray4D<real> Flux = data->hydro.FluxRiemann;
  IdefixArray1D<int> vars = this->varList;
  idefix_for("RKL_ResetFlux",
             0,nvarRKL,
             0,data->np_tot[KDIR],
             0,data->np_tot[JDIR],
             0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int n, int k, int j, int i) {
      const int nv = vars(n);
      Flux(nv,k,j,i) = ZERO_F;
    }
  );
  idfx::popRegion();
}

void RKLegendre::ResetStage() {
  idfx::pushRegion("RKLegendre::ResetStage");

  IdefixArray4D<real> dU = this->dU;
  IdefixArray4D<real> Flux = data->hydro.FluxRiemann;
  IdefixArray4D<real> dB = this->dB;
  IdefixArray3D<real> ex = data->hydro.emf.ex;
  IdefixArray3D<real> ey = data->hydro.emf.ey;
  IdefixArray3D<real> ez = data->hydro.emf.ez;
  IdefixArray1D<int> vars = this->varList;
  IdefixArray3D<real> invDt = data->hydro.InvDt;
  int stage = this->stage;
  int nvar = this->nvarRKL;
  bool haveVs=this->haveVs;

  idefix_for("RKL_ResetStage",
             0,data->np_tot[KDIR],
             0,data->np_tot[JDIR],
             0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      for(int n = 0 ; n < nvar ; n++) {
        const int nv = vars(n);
        dU(nv,k,j,i) = ZERO_F;
      }
      if(stage == 1)
        invDt(k,j,i) = ZERO_F;

      if(haveVs) {
        for(int n=0; n < DIMENSIONS; n++) {
          dB(n,k,j,i) = ZERO_F;
        }
        D_EXPAND( ez(k,j,i) = 0.0;    ,
                                      ,
                  ex(k,j,i) = 0.0;
                  ey(k,j,i) = 0.0;    )
      }
    });

  idfx::popRegion();
}


void RKLegendre::ComputeDt() {
  idfx::pushRegion("RKLegendre::ComputeDt");

  IdefixArray3D<real> invDt = data->hydro.InvDt;

  real newinvdt = ZERO_F;
  Kokkos::parallel_reduce("RKL_Timestep_reduction",
    Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>
    ({0,0,0},{data->end[KDIR],data->end[JDIR],data->end[IDIR]}),
    KOKKOS_LAMBDA (int k, int j, int i, real &invdt) {
      invdt = std::fmax(invDt(k,j,i), invdt);
    },
    Kokkos::Max<real>(newinvdt)
  );

#ifdef WITH_MPI
  if(idfx::psize>1) {
    MPI_SAFE_CALL(MPI_Allreduce(MPI_IN_PLACE, &newinvdt, 1, realMPI, MPI_MAX, MPI_COMM_WORLD));
  }
#endif

  dt = 1.0/newinvdt;
  dt = (cfl_rkl*dt)/2.0; // parabolic time step

  idfx::popRegion();
}

template<int dir> void RKLegendre::LoopDir(real t) {
    ResetFlux();

    // CalcParabolicFlux
    data->hydro.CalcParabolicFlux<dir>(t);

    // Calc Right Hand Side
    CalcParabolicRHS<dir>(t);

    // Recursive: do next dimension
    LoopDir<dir+1>(t);
}

template<> void RKLegendre::LoopDir<DIMENSIONS>(real t) {
  // Do nothing
}

void RKLegendre::EvolveStage(real t) {
  idfx::pushRegion("RKLegendre::EvolveStage");

  ResetStage();

  if(haveVs && data->hydro.needRKLCurrent) data->hydro.CalcCurrent();

  // Loop on dimensions for the parabolic fluxes and RHS, starting from IDIR
  LoopDir<IDIR>(t);

  if(haveVs) {
    data->hydro.emf.CalcNonidealEMF(t);
    data->hydro.emf.EnforceEMFBoundary();
    real dt=1.0;
    data->hydro.emf.EvolveMagField(t, dt, this->dB);
  }
  idfx::popRegion();
}

template <int dir>
void RKLegendre::CalcParabolicRHS(real t) {
  idfx::pushRegion("RKLegendre::CalcParabolicRHS");

  IdefixArray4D<real> Flux = data->hydro.FluxRiemann;
  IdefixArray3D<real> A    = data->A[dir];
  IdefixArray3D<real> dV   = data->dV;
  IdefixArray1D<real> x1m  = data->xl[IDIR];
  IdefixArray1D<real> x1   = data->x[IDIR];
  IdefixArray1D<real> sm   = data->sinx2m;
  IdefixArray1D<real> rt   = data->rt;
  IdefixArray1D<real> dmu  = data->dmu;
  IdefixArray1D<real> s    = data->sinx2;
  IdefixArray1D<real> dx   = data->dx[dir];
  IdefixArray1D<real> dx2  = data->dx[JDIR];
  IdefixArray3D<real> invDt = data->hydro.InvDt;
  IdefixArray3D<real> dMax = data->hydro.dMax;
  IdefixArray4D<real> viscSrc = data->hydro.viscosity.viscSrc;
  IdefixArray4D<real> dU = this->dU;
  IdefixArray1D<int> varList = this->varList;

  bool haveViscosity = data->hydro.viscosityStatus.isRKL;

  int ioffset,joffset,koffset;
  ioffset=joffset=koffset=0;
  // Determine the offset along which we do the extrapolation
  if(dir==IDIR) ioffset=1;
  if(dir==JDIR) joffset=1;
  if(dir==KDIR) koffset=1;


  idefix_for("CalcTotalFlux",
             0, this->nvarRKL,
             data->beg[KDIR],data->end[KDIR]+koffset,
             data->beg[JDIR],data->end[JDIR]+joffset,
             data->beg[IDIR],data->end[IDIR]+ioffset,
    KOKKOS_LAMBDA (int n, int k, int j, int i) {
      real Ax = A(k,j,i);

#if GEOMETRY != CARTESIAN
      if(Ax<SMALL_NUMBER)
        Ax=SMALL_NUMBER;    // Essentially to avoid singularity around poles
#endif

      const int nv = varList(n);

      Flux(nv,k,j,i) = Flux(nv,k,j,i) * Ax;

      // Curvature terms
#if    (GEOMETRY == POLAR       && COMPONENTS >= 2) \
    || (GEOMETRY == CYLINDRICAL && COMPONENTS == 3)
      if(dir==IDIR && nv==iMPHI) {
        // Conserve angular momentum, hence flux is R*Vphi
        Flux(iMPHI,k,j,i) = Flux(iMPHI,k,j,i) * FABS(x1m(i));
      }
#endif // GEOMETRY==POLAR OR CYLINDRICAL

#if GEOMETRY == SPHERICAL
      if(dir==IDIR && nv==iMPHI) {
  #if COMPONENTS == 3
        Flux(iMPHI,k,j,i) = Flux(iMPHI,k,j,i) * FABS(x1m(i));
  #endif // COMPONENTS == 3
      } else if(dir==JDIR && nv==iMPHI) {
  #if COMPONENTS == 3
        Flux(iMPHI,k,j,i) = Flux(iMPHI,k,j,i) * FABS(sm(j));
  #endif // COMPONENTS = 3
      }
#endif // GEOMETRY == SPHERICAL
    }
  );


  idefix_for("CalcRightHandSide",
             0, this->nvarRKL,
             data->beg[KDIR],data->end[KDIR],
             data->beg[JDIR],data->end[JDIR],
             data->beg[IDIR],data->end[IDIR],
    KOKKOS_LAMBDA (int n, int k, int j, int i) {
      constexpr const int ioffset = (dir==IDIR) ? 1 : 0;
      constexpr const int joffset = (dir==JDIR) ? 1 : 0;
      constexpr const int koffset = (dir==KDIR) ? 1 : 0;
      real rhs;

      const int nv = varList(n);

      rhs = -  ( Flux(nv, k+koffset, j+joffset, i+ioffset)
                     - Flux(nv, k, j, i))/dV(k,j,i);

      // Viscosity source terms
      if( haveViscosity && (nv-VX1 < COMPONENTS) && (nv-VX1>=0)) {
        rhs += viscSrc(nv-VX1,k,j,i);
      }

#if GEOMETRY != CARTESIAN
  #ifdef iMPHI
      if((dir==IDIR) && (nv == iMPHI)) {
        rhs /= x1(i);
      }
    #if (GEOMETRY == SPHERICAL) && (COMPONENTS == 3)
      if((dir==JDIR) && (nv == iMPHI)) {
        rhs /= FABS(s(j));
      }
    #endif // GEOMETRY
      // Nothing for KDIR
  #endif  // iMPHI
#endif // GEOMETRY != CARTESIAN

      // store the field components
      dU(nv,k,j,i) += rhs;
    });

  // Compute hyperbolic timestep only if we're in the first stage of the RKL loop
  if(stage==1) {
    idefix_for("CalcDt",
             data->beg[KDIR],data->end[KDIR],
             data->beg[JDIR],data->end[JDIR],
             data->beg[IDIR],data->end[IDIR],
             KOKKOS_LAMBDA (int k, int j, int i) {
                constexpr const int ioffset = (dir==IDIR) ? 1 : 0;
                constexpr const int joffset = (dir==JDIR) ? 1 : 0;
                constexpr const int koffset = (dir==KDIR) ? 1 : 0;
               // Compute dt from max signal speed
                const int ig = ioffset*i + joffset*j + koffset*k;
                real dl = dx(ig);
                #if GEOMETRY == POLAR
                  if (dir==JDIR)
                    dl = dl*x1(i);

                #elif GEOMETRY == SPHERICAL
                  if (dir==JDIR)
                    dl = dl*rt(i);
                  else
                    if (dir==KDIR)
                      dl = dl*rt(i)*dmu(j)/dx2(j);
                 #endif

                invDt(k,j,i) += 0.5 * std::fmax(dMax(k+koffset,j+joffset,i+ioffset),
                                                        dMax(k,j,i)) / (dl*dl);
      });
  }

  idfx::popRegion();
}

void RKLegendre::SetBoundaries(real t) {
  idfx::pushRegion("RKLegendre::SetBoundaries");
  // set internal boundary conditions
  if(data->hydro.boundary.haveInternalBoundary) data->hydro.boundary.internalBoundaryFunc(*data, t);
  for(int dir=0 ; dir < DIMENSIONS ; dir++ ) {
      // MPI Exchange data when needed
      // We use the RKL instance MPI object to ensure that we only exchange the data
      // solved by RKL
    #ifdef WITH_MPI
    if(data->mygrid->nproc[dir]>1) {
      switch(dir) {
        case 0:
          this->mpi.ExchangeX1();
          break;
        case 1:
          this->mpi.ExchangeX2();
          break;
        case 2:
          this->mpi.ExchangeX3();
          break;
      }
    }
    #endif
    data->hydro.boundary.EnforceBoundaryDir(t, dir);
    #if MHD == YES
      // Reconstruct the normal field component when using CT
      if(haveVs) {
        data->hydro.boundary.ReconstructNormalField(dir);
      }
    #endif
  } // Loop on dimension ends

#if MHD == YES
  // Remake the cell-centered field.
  if(haveVs) {
    data->hydro.boundary.ReconstructVcField(data->hydro.Vc);
  }
#endif
  idfx::popRegion();
}
