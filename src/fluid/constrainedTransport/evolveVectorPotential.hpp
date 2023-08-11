// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_CONSTRAINEDTRANSPORT_EVOLVEVECTORPOTENTIAL_HPP_
#define FLUID_CONSTRAINEDTRANSPORT_EVOLVEVECTORPOTENTIAL_HPP_

#include <string>

#include "idefix.hpp"
#include "fluid.hpp"
#include "dataBlock.hpp"

template<typename Phys>
void ConstrainedTransport<Phys>::EvolveVectorPotential(real dt, IdefixArray4D<real> &Vein) {
  #ifdef EVOLVE_VECTOR_POTENTIAL
    idfx::pushRegion("ConstrainedTransport::EvolveVectorPotential");
    IdefixArray4D<real> Ve = Vein;
        // Corned EMFs
    IdefixArray3D<real> Ex1 = this->ex;
    IdefixArray3D<real> Ex2 = this->ey;
    IdefixArray3D<real> Ex3 = this->ez;
    idefix_for("EvolvVectorPotential",
      data->beg[KDIR],data->end[KDIR]+KOFFSET,
      data->beg[JDIR],data->end[JDIR]+JOFFSET,
      data->beg[IDIR],data->end[IDIR]+IOFFSET,
      KOKKOS_LAMBDA (int k, int j, int i) {
        #if DIMENSIONS == 3
          Ve(AX1e,k,j,i) += - dt * Ex1(k,j,i);
          Ve(AX2e,k,j,i) += - dt * Ex2(k,j,i);
        #endif
        #if DIMENSIONS >= 2
          Ve(AX3e,k,j,i) += - dt * Ex3(k,j,i);
        #endif
      });

    idfx::popRegion();
  #endif
}



template<typename Phys>
void ConstrainedTransport<Phys>::ComputeMagFieldFromA(IdefixArray4D<real> &Vein,
                                              IdefixArray4D<real> &Vsout) {
  #ifdef EVOLVE_VECTOR_POTENTIAL
    idfx::pushRegion("ConstrainedTransport::ComputeMagFieldfromA");

    // Corned EMFs
    IdefixArray3D<real> Ex1 = this->ex;
    IdefixArray3D<real> Ex2 = this->ey;
    IdefixArray3D<real> Ex3 = this->ez;

    // Field
    IdefixArray4D<real> Vs = Vsout;
    IdefixArray4D<real> Ve = Vein;

    // Coordinates
    IdefixArray1D<real> x1=data->x[IDIR];
    IdefixArray1D<real> x2=data->x[JDIR];
    IdefixArray1D<real> x3=data->x[KDIR];

    IdefixArray1D<real> x1p=data->xr[IDIR];
    IdefixArray1D<real> x2p=data->xr[JDIR];
    IdefixArray1D<real> x3p=data->xr[KDIR];

    IdefixArray1D<real> x1m=data->xl[IDIR];
    IdefixArray1D<real> x2m=data->xl[JDIR];
    IdefixArray1D<real> x3m=data->xl[KDIR];

    IdefixArray1D<real> dx1=data->dx[IDIR];
    IdefixArray1D<real> dx2=data->dx[JDIR];
    IdefixArray1D<real> dx3=data->dx[KDIR];

    #if GEOMETRY == SPHERICAL
      IdefixArray1D<real> dmu=data->dmu;
      IdefixArray1D<real> sinx2m=data->sinx2m;
      #if DIMENSIONS >= 2
        bool haveAxis = hydro->haveAxis;
      #endif
    #endif



    idefix_for("ComputeMagFieldFromA",
              data->beg[KDIR],data->end[KDIR]+KOFFSET,
              data->beg[JDIR],data->end[JDIR]+JOFFSET,
              data->beg[IDIR],data->end[IDIR]+IOFFSET,
      KOKKOS_LAMBDA (int k, int j, int i) {
  #if GEOMETRY == CARTESIAN
    Vs(BX1s,k,j,i) = D_EXPAND( ZERO_F                                     ,
                    + 1/dx2(j) * (Ve(AX3e,k,j+1,i) - Ve(AX3e,k,j,i) )  ,
                    - 1/dx3(k) * (Ve(AX2e,k+1,j,i) - Ve(AX2e,k,j,i) )  );

    #if DIMENSIONS >= 2
      Vs(BX2s,k,j,i) =  D_EXPAND( - 1/dx1(i) * (Ve(AX3e,k,j,i+1) - Ve(AX3e,k,j,i) )  ,
                                                                  ,
                        + 1/dx3(k) * (Ve(AX1e,k+1,j,i) - Ve(AX1e,k,j,i) ) );
    #endif
    #if DIMENSIONS == 3
      Vs(BX3s,k,j,i) = + 1/dx1(i) * (Ve(AX2e,k,j,i+1) - Ve(AX2e,k,j,i) )
              - 1/dx2(j) * (Ve(AX1e,k,j+1,i) - Ve(AX1e,k,j,i) );
    #endif

  #elif GEOMETRY == CYLINDRICAL
    Vs(BX1s,k,j,i) = + 1/dx2(j) * (Ve(AX3e,k,j+1,i) - Ve(AX3e,k,j,i) );
    #if DIMENSIONS >= 2
      Vs(BX2s,k,j,i) =  -( FABS(x1p(i)) * Ve(AX3e,k,j,i+1)
                          - FABS(x1m(i)) * Ve(AX3e,k,j,i)) / FABS(x1(i)*dx1(i));
    #endif

  #elif GEOMETRY == POLAR
    Vs(BX1s,k,j,i) = D_EXPAND( ZERO_F                                                      ,
                    + 1/(FABS(x1m(i)) * dx2(j)) * (Ve(AX3e,k,j+1,i) - Ve(AX3e,k,j,i) )  ,
                    - 1/dx3(k) * (Ve(AX2e,k+1,j,i) - Ve(AX2e,k,j,i) )                   );

    #if DIMENSIONS >= 2
      Vs(BX2s,k,j,i) =  D_EXPAND( - 1/dx1(i) * (Ve(AX3e,k,j,i+1) - Ve(AX3e,k,j,i) )  ,
                                                                  ,
                        + 1/dx3(k) * (Ve(AX1e,k+1,j,i) - Ve(AX1e,k,j,i) ) );
    #endif
    #if DIMENSIONS == 3
      Vs(BX3s,k,j,i) = 1/(FABS(x1(i))) * (
                    (x1m(i+1)*Ve(AX2e,k,j,i+1) - x1m(i)*Ve(AX2e,k,j,i) ) / dx1(i)
                  -  (Ve(AX1e,k,j+1,i) - Ve(AX1e,k,j,i) ) / dx2(j) );
    #endif
  #elif GEOMETRY == SPHERICAL
    real dV2  = dmu(j);
    real Ax2p = FABS(sinx2m(j+1));
    real Ax2m = FABS(sinx2m(j));

    Vs(BX1s,k,j,i) = D_EXPAND( ZERO_F                                                        ,
                    + 1/(x1m(i)*dV2) * ( Ax2p*Ve(AX3e,k,j+1,i) - Ax2m*Ve(AX3e,k,j,i) )    ,
                    - dx2(j)/(x1m(i)*dV2*dx3(k)) * (Ve(AX2e,k+1,j,i) - Ve(AX2e,k,j,i) ) );

    #if DIMENSIONS >= 2
      // If we include the axis, we symmetrize Ex on the axis. However, Ax2=0 on the axis
      // so BX2s might become singular. We therefore enforce Ax2=1 on the axis, knowing
      // that the contribution to BX2s of this term will be zero because Ve(AX1e,k+1)-Ve(AX1e,k)=0
      if(haveAxis) {
        if(FABS(Ax2m)<1e-12) Ax2m = ONE_F;
      }
      Vs(BX2s,k,j,i) =  D_EXPAND( - 1/(x1(i)*dx1(i)) * (x1m(i+1)*Ve(AX3e,k,j,i+1)
                                                        - x1m(i)*Ve(AX3e,k,j,i) )  ,
                                                                                  ,
                        + 1/(x1(i)*Ax2m*dx3(k)) * (Ve(AX1e,k+1,j,i) - Ve(AX1e,k,j,i) )  );
    #endif
    #if DIMENSIONS == 3
      Vs(BX3s,k,j,i) = + 1/(x1(i)*dx1(i)) * (x1m(i+1)*Ve(AX2e,k,j,i+1) - x1m(i)*Ve(AX2e,k,j,i) )
              - 1/(x1(i)*dx2(j)) * (Ve(AX1e,k,j+1,i) - Ve(AX1e,k,j,i) );
    #endif
  #endif // GEOMETRY
  });

  if(data->haveGridCoarsening) {
    hydro->CoarsenMagField(Vsout);
  }

  idfx::popRegion();
  #endif
}
#endif //FLUID_CONSTRAINEDTRANSPORT_EVOLVEVECTORPOTENTIAL_HPP_
