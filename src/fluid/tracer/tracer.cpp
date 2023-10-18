// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "idefix.hpp"
#include "tracer.hpp"

// Convect Conservative to Primitive variable
void Tracer::ConvertConsToPrim() {
  idfx::pushRegion("Tracer::ConvertConsToPrim");

  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Uc = this->Uc;

  idefix_for("ConsToPrimScalar",
            nVar, nVar+nTracer,   // Loop on the index where scalars are lying
             0,data->np_tot[KDIR],
             0,data->np_tot[JDIR],
             0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int n, int k, int j, int i) {
      Vc(n,k,j,i) = Uc(n,k,j,i) / Uc(RHO,k,j,i);
  });

  idfx::popRegion();
}



// Convect Conservative to Primitive variable
void Tracer::ConvertPrimToCons() {
  idfx::pushRegion("Tracer::ConvertPrimToCons");

  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Uc = this->Uc;

  idefix_for("PrimToConsScalar",
             nVar, nVar+nTracer,  // Loop on the index where scalars are lying
             0,data->np_tot[KDIR],
             0,data->np_tot[JDIR],
             0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int n, int k, int j, int i) {
      Uc(n,k,j,i) = Vc(n,k,j,i) * Vc(RHO,k,j,i);
  });

  idfx::popRegion();
}
