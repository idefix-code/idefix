// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_CONSTRAINEDTRANSPORT_CALCNONIDEALEMF_HPP_
#define FLUID_CONSTRAINEDTRANSPORT_CALCNONIDEALEMF_HPP_

#include "fluid.hpp"
#include "dataBlock.hpp"

// Compute Corner EMFs from nonideal MHD
template<typename Phys>
void ConstrainedTransport<Phys>::CalcNonidealEMF(real t) {
  idfx::pushRegion("ConstrainedTransport::CalcNonidealEMF");

#if MHD == YES
  // Corned EMFs
  IdefixArray3D<real> ex = this->ex;
  IdefixArray3D<real> ey = this->ey;
  IdefixArray3D<real> ez = this->ez;
  IdefixArray4D<real> J = hydro->J;
  IdefixArray4D<real> Vs = hydro->Vs;
  IdefixArray4D<real> Vc = hydro->Vc;

  // These arrays have been previously computed in calcParabolicFlux
  IdefixArray3D<real> etaArr = hydro->etaOhmic;
  IdefixArray3D<real> xAmbiArr = hydro->xAmbipolar;

  // these two are required to ensure that the type is captured by KOKKOS_LAMBDA
  HydroModuleStatus resistivity = hydro->resistivityStatus.status;
  HydroModuleStatus ambipolar = hydro->ambipolarStatus.status;

  bool haveResistivity{false};
  bool haveAmbipolar{false};

  if(data->rklCycle) {
    haveResistivity = hydro->resistivityStatus.isRKL;
    haveAmbipolar = hydro->ambipolarStatus.isRKL;
  } else {
    haveResistivity = hydro->resistivityStatus.isExplicit;
    haveAmbipolar = hydro->ambipolarStatus.isExplicit;
  }

  real etaConstant = hydro->etaO;
  real xAConstant = hydro->xA;

  idefix_for("CalcNIEMF",
             data->beg[KDIR],data->end[KDIR]+KOFFSET,
             data->beg[JDIR],data->end[JDIR]+JOFFSET,
             data->beg[IDIR],data->end[IDIR]+IOFFSET,
    KOKKOS_LAMBDA (int k, int j, int i) {
      real Bx1, Bx2, Bx3;
      real Jx1, Jx2, Jx3;
      real eta, xA;
      // CT_EMF_ArithmeticAverage (emf, 0.25);

      if(resistivity == Constant)
        eta = etaConstant;
      if(ambipolar == Constant)
        xA = xAConstant;

  #if DIMENSIONS == 3
      // -----------------------
      // X1 EMF Component
      // -----------------------
      Jx1 = J(IDIR,k,j,i);

      // Ohmic resistivity

      if(haveResistivity) {
        if(resistivity == UserDefFunction) eta = AVERAGE_3D_YZ(etaArr,k,j,i);
        ex(k,j,i) += eta * Jx1;
      }

      // Ambipolar diffusion
      if(haveAmbipolar) {
        if(ambipolar == UserDefFunction) xA = AVERAGE_3D_YZ(xAmbiArr,k,j,i);
        Bx1 = AVERAGE_4D_XYZ(Vs, BX1s, k,j,i+1);
        Bx2 = AVERAGE_4D_Z(Vs, BX2s, k, j, i);
        Bx3 = AVERAGE_4D_Y(Vs, BX3s, k, j, i);

        // Jx1 is already defined above
        Jx2 = AVERAGE_4D_XY(J, JDIR, k, j, i+1);
        Jx3 = AVERAGE_4D_XZ(J, KDIR, k, j, i+1);

        real JdotB = (Jx1*Bx1 + Jx2*Bx2 + Jx3*Bx3);
        real BdotB = (Bx1*Bx1 + Bx2*Bx2 + Bx3*Bx3);

        ex(k,j,i) += xA * (BdotB*Jx1 - JdotB * Bx1);
      }

      // -----------------------
      // X2 EMF Component
      // -----------------------
      Jx2 = J(JDIR,k,j,i);

      // Ohmic resistivity
      if(haveResistivity) {
        if(resistivity == UserDefFunction) eta = AVERAGE_3D_XZ(etaArr,k,j,i);
        ey(k,j,i) += eta * Jx2;
      }

      // Ambipolar diffusion
      if(haveAmbipolar) {
        if(ambipolar == UserDefFunction) xA = AVERAGE_3D_XZ(xAmbiArr,k,j,i);
        Bx1 = AVERAGE_4D_Z(Vs, BX1s, k, j, i);
        Bx2 = AVERAGE_4D_XYZ(Vs, BX2s, k, j+1, i);
        Bx3 = AVERAGE_4D_X(Vs, BX3s, k, j, i);

        // Jx2 is already defined above
        Jx1 = AVERAGE_4D_XY(J, IDIR, k, j+1, i);
        Jx3 = AVERAGE_4D_YZ(J, KDIR, k, j+1, i);

        real JdotB = (Jx1*Bx1 + Jx2*Bx2 + Jx3*Bx3);
        real BdotB = (Bx1*Bx1 + Bx2*Bx2 + Bx3*Bx3);

        ey(k,j,i) += xA * (BdotB*Jx2 - JdotB * Bx2);
      }
  #endif
      // -----------------------
      // X3 EMF Component
      // -----------------------
      Jx3 = J(KDIR,k,j,i);

      // Ohmic resistivity
      if(haveResistivity) {
        if(resistivity == UserDefFunction) eta = AVERAGE_3D_XY(etaArr,k,j,i);
        ez(k,j,i) += eta * Jx3;
      }

      // Ambipolar diffusion
      if(haveAmbipolar) {
        if(ambipolar == UserDefFunction) xA = AVERAGE_3D_XY(xAmbiArr,k,j,i);
        Bx1 = AVERAGE_4D_Y(Vs, BX1s, k, j, i);
  #if DIMENSIONS >= 2
        Bx2 = AVERAGE_4D_X(Vs, BX2s, k, j, i);
  #else
    #if COMPONENTS >= 2
        Bx2 = AVERAGE_4D_XY(Vc, BX2, k, j, i);
    #else
        Bx2 = 0.0;
    #endif
  #endif

  #if DIMENSIONS == 3
        Bx3 = AVERAGE_4D_XYZ(Vs, BX3s, k+1, j, i);
  #else
      #if COMPONENTS == 3
        Bx3 = AVERAGE_4D_XY(Vc, BX3, k, j, i);
      #else
        Bx3 = 0.0;
      #endif
  #endif

        // Jx3 is already defined above
  #if DIMENSIONS == 3
        Jx1 = AVERAGE_4D_XZ(J, IDIR, k+1, j, i);
        Jx2 = AVERAGE_4D_YZ(J, JDIR, k+1, j, i);
  #else
        Jx1 = AVERAGE_4D_X(J, IDIR, k, j, i);
        Jx2 = AVERAGE_4D_Y(J, JDIR, k, j, i);
  #endif
        real JdotB = (Jx1*Bx1 + Jx2*Bx2 + Jx3*Bx3);
        real BdotB = (Bx1*Bx1 + Bx2*Bx2 + Bx3*Bx3);

        ez(k,j,i) += xA * (BdotB * Jx3 - JdotB * Bx3);
      }
    }
  );
#endif

  idfx::popRegion();
}

#endif //FLUID_CONSTRAINEDTRANSPORT_CALCNONIDEALEMF_HPP_
