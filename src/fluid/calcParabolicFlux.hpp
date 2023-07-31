// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************
#ifndef FLUID_CALCPARABOLICFLUX_HPP_
#define FLUID_CALCPARABOLICFLUX_HPP_

#include "fluid.hpp"
#include "dataBlock.hpp"
#include "addNonIdealMHDFlux.hpp"
#include "fargo.hpp"

// Compute parabolic fluxes
template <typename Phys>
template <int dir>
void Fluid<Phys>::CalcParabolicFlux(const real t) {
  idfx::pushRegion("Fluid::CalcParabolicFlux");

  IdefixArray3D<real> dMax = this->dMax;

  // Reset Max diffusion coefficient
  idefix_for("HydroParabolicResetStage",0,data->np_tot[KDIR],
                                        0,data->np_tot[JDIR],
                                        0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
        dMax(k,j,i) = ZERO_F;
      }
  );

  if( (resistivityStatus.isExplicit  && (! data->rklCycle))
    || (resistivityStatus.isRKL  && ( data->rklCycle))
    || (ambipolarStatus.isExplicit  && (! data->rklCycle))
    || (ambipolarStatus.isRKL  && ( data->rklCycle)) ) {
      this->AddNonIdealMHDFlux<dir>(t);
  }

  if( (viscosityStatus.isExplicit && (!data->rklCycle))
    || (viscosityStatus.isRKL && data->rklCycle))  {
      // Add fargo velocity if using fargo
    if(data->haveFargo && viscosityStatus.isExplicit) {
      data->fargo->AddVelocityFluid(t,this);
    }
    this->viscosity->AddViscousFlux(dir,t, this->FluxRiemann);

    // Remove back Fargo velocity
    if(data->haveFargo && viscosityStatus.isExplicit) {
      data->fargo->SubstractVelocityFluid(t,this);
    }
  }

  // Add thermal diffusion
  if( (thermalDiffusionStatus.isExplicit && (!data->rklCycle))
    || (thermalDiffusionStatus.isRKL && data->rklCycle))  {
    this->thermalDiffusion->AddDiffusiveFlux(dir,t, this->FluxRiemann);
  }

  if( (bragViscosityStatus.isExplicit && (!data->rklCycle))
    || (bragViscosityStatus.isRKL && data->rklCycle))  {
    this->bragViscosity->AddBragViscousFlux(dir,t, this->FluxRiemann);
  }

  // Add braginskii thermal diffusion
  if( (bragThermalDiffusionStatus.isExplicit && (!data->rklCycle))
    || (bragThermalDiffusionStatus.isRKL && data->rklCycle))  {
    this->bragThermalDiffusion->AddBragDiffusiveFlux(dir,t, this->FluxRiemann);
  }

  idfx::popRegion();
}

#endif //FLUID_CALCPARABOLICFLUX_HPP_
