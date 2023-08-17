// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_CONVERTCONSTOPRIM_HPP_
#define FLUID_CONVERTCONSTOPRIM_HPP_

#include "fluid.hpp"
#include "dataBlock.hpp"
#include "tracer.hpp"

template <typename Phys>
KOKKOS_INLINE_FUNCTION void K_ConsToPrim(real Vc[], real Uc[], const EquationOfState *eos) {
  Vc[RHO] = Uc[RHO];

  EXPAND( Vc[VX1] = Uc[MX1]/Uc[RHO];  ,
          Vc[VX2] = Uc[MX2]/Uc[RHO];  ,
          Vc[VX3] = Uc[MX3]/Uc[RHO];  )

  if constexpr(Phys::mhd) {
    EXPAND( Vc[BX1] = Uc[BX1];  ,
            Vc[BX2] = Uc[BX2];  ,
            Vc[BX3] = Uc[BX3];  )
  }


  if constexpr(Phys::pressure) {
    real kin = HALF_F / Uc[RHO] * (EXPAND( Uc[MX1]*Uc[MX1]   ,
                                    + Uc[MX2]*Uc[MX2]  ,
                                    + Uc[MX3]*Uc[MX3]  ));

    if constexpr(Phys::mhd) {
      real mag = HALF_F * (EXPAND( Uc[BX1]*Uc[BX1]   ,
                          + Uc[BX2]*Uc[BX2]  ,
                          + Uc[BX3]*Uc[BX3]  ));

      Vc[PRS] = eos->GetPressure(Uc[ENG] - kin - mag, Uc[RHO]);

      // Check pressure positivity
      if(Vc[PRS]<= ZERO_F) {
        #ifdef SMALL_PRESSURE_TEMPERATURE
          Vc[PRS] = SMALL_PRESSURE_TEMPERATURE*Vc[RHO];
        #else
          Vc[PRS] = SMALL_PRESSURE_FIX;
        #endif

          Uc[ENG] = eos->GetInternalEnergy(Vc[PRS],Vc[RHO]) + kin + mag;
      }

    } else { // Hydro case
      Vc[PRS] = eos->GetPressure(Uc[ENG] - kin, Uc[RHO]);
      // Check pressure positivity
      if(Vc[PRS]<= ZERO_F) {
        #ifdef SMALL_PRESSURE_TEMPERATURE
          Vc[PRS] = SMALL_PRESSURE_TEMPERATURE*Vc[RHO];
        #else
          Vc[PRS] = SMALL_PRESSURE_FIX;
        #endif

          Uc[ENG] = eos->GetInternalEnergy(Vc[PRS],Vc[RHO]) + kin;
      }
    } // MHD
  } // Have Energy
}

template <typename Phys>
KOKKOS_INLINE_FUNCTION void K_PrimToCons(real Uc[], real Vc[], const EquationOfState *eos) {
  Uc[RHO] = Vc[RHO];

  EXPAND( Uc[MX1] = Vc[VX1]*Vc[RHO];  ,
          Uc[MX2] = Vc[VX2]*Vc[RHO];  ,
          Uc[MX3] = Vc[VX3]*Vc[RHO];  )

  if constexpr(Phys::mhd) {
    EXPAND( Uc[BX1] = Vc[BX1];  ,
            Uc[BX2] = Vc[BX2];  ,
            Uc[BX3] = Vc[BX3];  )
  }


  if constexpr(Phys::pressure) {
    if constexpr(Phys::mhd) {
      Uc[ENG] = eos->GetInternalEnergy(Vc[PRS],Vc[RHO])
                + HALF_F * Vc[RHO] * (EXPAND( Vc[VX1]*Vc[VX1]  ,
                                            + Vc[VX2]*Vc[VX2]  ,
                                            + Vc[VX3]*Vc[VX3]  ))
                + HALF_F * (EXPAND( Uc[BX1]*Uc[BX1]  ,
                                  + Uc[BX2]*Uc[BX2]  ,
                                  + Uc[BX3]*Uc[BX3]  ));
    } else {
      Uc[ENG] = eos->GetInternalEnergy(Vc[PRS],Vc[RHO])
                + HALF_F * Vc[RHO] * (EXPAND( Vc[VX1]*Vc[VX1]  ,
                                            + Vc[VX2]*Vc[VX2]  ,
                                            + Vc[VX3]*Vc[VX3]  ));
    } //MHD
  } // Energy
}



// Convect Conservative to Primitive variable
template<typename Phys>
void Fluid<Phys>::ConvertConsToPrim() {
  idfx::pushRegion("Fluid::ConvertConsToPrim");

  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Uc = this->Uc;
  EquationOfState eos;
  if constexpr(Phys::eos) {
    eos = *(this->eos.get());
  }

  if constexpr(Phys::mhd) {
    #ifdef EVOLVE_VECTOR_POTENTIAL
      emf->ComputeMagFieldFromA(Ve,Vs);
    #endif
    boundary->ReconstructVcField(Uc);
  }

  idefix_for("ConsToPrim",
             0,data->np_tot[KDIR],
             0,data->np_tot[JDIR],
             0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      real U[Phys::nvar];
      real V[Phys::nvar];

#pragma unroll
      for(int nv = 0 ; nv < Phys::nvar; nv++) {
        U[nv] = Uc(nv,k,j,i);
      }

      K_ConsToPrim<Phys>(V,U,&eos);

#pragma unroll
      for(int nv = 0 ; nv<Phys::nvar; nv++) {
        Vc(nv,k,j,i) = V[nv];
      }
  });

  if(haveTracer) {
    tracer->ConvertConsToPrim();
  }

  idfx::popRegion();
}

// Convert Primitive to conservative variables
template<typename Phys>
void Fluid<Phys>::ConvertPrimToCons() {
  idfx::pushRegion("Fluid::ConvertPrimToCons");

  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Uc = this->Uc;
  EquationOfState eos;
  if constexpr(Phys::eos) {
    eos = *(this->eos.get());
  }

  idefix_for("ConvertPrimToCons",
             0,data->np_tot[KDIR],
             0,data->np_tot[JDIR],
             0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      real U[Phys::nvar];
      real V[Phys::nvar];

#pragma unroll
      for(int nv = 0 ; nv < Phys::nvar; nv++) {
        V[nv] = Vc(nv,k,j,i);
      }

      K_PrimToCons<Phys>(U,V,&eos);

#pragma unroll
      for(int nv = 0 ; nv<Phys::nvar; nv++) {
        Uc(nv,k,j,i) = U[nv];
      }
  });

  if(haveTracer) {
    tracer->ConvertPrimToCons();
  }

  idfx::popRegion();
}

#endif //FLUID_CONVERTCONSTOPRIM_HPP_
