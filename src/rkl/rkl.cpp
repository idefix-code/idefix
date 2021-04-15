// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <cmath>

#include "rkl.hpp"
#include "rkl_defs.hpp"
#include "dataBlock.hpp"


RKLegendre::RKLegendre() {
  // do nothing!
}


void RKLegendre::Init(DataBlock *datain) {
  idfx::pushRegion("RKLegendre::Init");

  // Save the datablock to which we are attached from now on
  this->data = datain;

  dU = IdefixArray4D<real>("RKL_dU", NVAR,
                           data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  dU0 = IdefixArray4D<real>("RKL_dU0", NVAR,
                           data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  Uc1 = IdefixArray4D<real>("RKL_Uc1", NVAR,
                           data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

  cfl_par = HALF_F;    //(no need for dimension since the dt definition in idefix is different)
  rmax_par = 100.0;

  idfx::popRegion();
}


void RKLegendre::Cycle() {
  idfx::pushRegion("RKLegendre::Cycle");

  real time = data->t;

  EvolveStage(time);

  idfx::popRegion();
}


void RKLegendre::ResetStage() {
  idfx::pushRegion("RKLegendre::ResetStage");

  IdefixArray4D<real> dU = this->dU;
  IdefixArray4D<real> Flux = data->hydro.FluxRiemann;
  idefix_for("InitRKLStage_dU_Flux",
             0,NVAR,
             0,data->np_tot[KDIR],
             0,data->np_tot[JDIR],
             0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int nv, int k, int j, int i) {
      dU(nv,k,j,i) = ZERO_F;
      Flux(nv,k,j,i) = ZERO_F;
    }
  );

  IdefixArray3D<real> invDt = data->hydro.InvDt;
  idefix_for("InitRKLStage_invDt",
             0,data->np_tot[KDIR],
             0,data->np_tot[JDIR],
             0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      invDt(k,j,i) = ZERO_F;
    }
  );

  this->dt = ZERO_F;

  idfx::popRegion();
}


void RKLegendre::ComputeDt() {
  idfx::pushRegion("RKLegendre::ComputeInvDt");

  IdefixArray3D<real> invDt = data->hydro.InvDt;

  real newinvdt = ZERO_F;
  Kokkos::parallel_reduce("Timestep_reduction",
    Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>
    ({data->beg[KDIR],data->beg[JDIR],data->beg[IDIR]},
    {data->end[KDIR],data->end[JDIR],data->end[IDIR]}),
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

  dt = ONE_F/newinvdt;
  dt = (cfl_par*dt)/TWO_F; // parabolic time step

  idfx::popRegion();
}


void RKLegendre::EvolveStage(real t) {
  idfx::pushRegion("RKLegendre::EvolveStage");

  ResetStage();

  for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
    // CalcParabolicFlux
    data->hydro.CalcParabolicFlux(dir, t);

    // Calc Right Hand Side
    CalcParabolicRHS(dir, t);
  }

  if (stage == 1) {
    ComputeDt();
  }

  idfx::popRegion();
}


void RKLegendre::CalcParabolicRHS(int dir, real t) {
  idfx::pushRegion("RKLegendre::CalcParabolicRHS");

  idfx::popRegion();
}
