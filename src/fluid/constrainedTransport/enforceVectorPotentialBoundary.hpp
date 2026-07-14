// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_CONSTRAINEDTRANSPORT_ENFORCEVECTORPOTENTIALBOUNDARY_HPP_
#define FLUID_CONSTRAINEDTRANSPORT_ENFORCEVECTORPOTENTIALBOUNDARY_HPP_
#include "constrainedTransport.hpp"

template<typename Phys>
void ConstrainedTransport<Phys>::EnforceVectorPotentialBoundary(IdefixArray4D<real> &Vein) {
  idfx::pushRegion("Emf::EnforceVectorPotentialBoundary");


  IdefixArray3D<real> Ax1, Ax2, Ax3;

  #ifdef EVOLVE_VECTOR_POTENTIAL
    #if DIMENSIONS == 3
      Ax1 = Kokkos::subview(Vein, AX1e, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
      Ax2 = Kokkos::subview(Vein, AX2e, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
    #endif
    Ax3 = Kokkos::subview(Vein, AX3e, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());


    if(this->hydro->haveAxis) {
      this->hydro->boundary->axis->RegularizeEMFs(Ax1, Ax2, Ax3);
    }

    #ifdef ENFORCE_EMF_CONSISTENCY
      #ifdef WITH_MPI
        // This average the vector potential at the domain surface with immediate neighbours
        // to ensure the vector potentials exactly match

        this->ExchangeAll(Ax1, Ax2, Ax3);
      #endif
      EnforceEMFBoundaryPeriodic(Ax1, Ax2, Ax3);
    #endif
  #endif // EVOLVE_VECTOR_POTENTIAL

  idfx::popRegion();
}
#endif // FLUID_CONSTRAINEDTRANSPORT_ENFORCEVECTORPOTENTIALBOUNDARY_HPP_
