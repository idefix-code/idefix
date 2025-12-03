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

  auto Ax1 = Kokkos::subview(Vein, IDIR, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
  auto Ax2 = Kokkos::subview(Vein, JDIR, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
  auto Ax3 = Kokkos::subview(Vein, KDIR, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

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

  idfx::popRegion();
}
#endif // FLUID_CONSTRAINEDTRANSPORT_ENFORCEVECTORPOTENTIALBOUNDARY_HPP_
