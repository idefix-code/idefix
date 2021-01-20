// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "hydro.hpp"
#include "dataBlock.hpp"

#if MHD == YES
#include "solversMHD.hpp"
#else
#include "solversHD.hpp"
#endif

// Convert Primitive to conservative variables
void Hydro::ConvertPrimToCons() {
  idfx::pushRegion("Hydro::ConvertPrimToCons");

  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Uc = this->Uc;
  real gamma_m1=this->gamma-ONE_F;

  idefix_for("ConvertPrimToCons",
             0,data->np_tot[KDIR],
             0,data->np_tot[JDIR],
             0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      real U[NVAR];
      real V[NVAR];

#pragma unroll
      for(int nv = 0 ; nv < NVAR; nv++) {
        V[nv] = Vc(nv,k,j,i);
      }

      K_PrimToCons(U,V,gamma_m1);

#pragma unroll
      for(int nv = 0 ; nv<NVAR; nv++) {
        Uc(nv,k,j,i) = U[nv];
      }
  });

  idfx::popRegion();
}
