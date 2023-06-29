// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_CHECKDIVB_HPP_
#define FLUID_CHECKDIVB_HPP_

#include "fluid.hpp"
#include "dataBlock.hpp"

template<typename Phys>
real Fluid<Phys>::CheckDivB() {
  real divB;
  IdefixArray4D<real> Vs = this->Vs;
  IdefixArray1D<real> dx1 = data->dx[IDIR];
  IdefixArray1D<real> dx2 = data->dx[JDIR];
  IdefixArray1D<real> dx3 = data->dx[KDIR];
  IdefixArray3D<real> Ax1 = data->A[IDIR];
  IdefixArray3D<real> Ax2 = data->A[JDIR];
  IdefixArray3D<real> Ax3 = data->A[KDIR];
  IdefixArray3D<real> dV = data->dV;


  idefix_reduce("CheckDivB",
    data->beg[KDIR], data->end[KDIR],
    data->beg[JDIR], data->end[JDIR],
    data->beg[IDIR], data->end[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i, real &divBmax) {
      [[maybe_unused]] real dB1,dB2,dB3;

      dB1=dB2=dB3=ZERO_F;

      D_EXPAND( dB1=(Ax1(k,j,i+1)*Vs(BX1s,k,j,i+1)-Ax1(k,j,i)*Vs(BX1s,k,j,i));  ,
                dB2=(Ax2(k,j+1,i)*Vs(BX2s,k,j+1,i)-Ax2(k,j,i)*Vs(BX2s,k,j,i));  ,
                dB3=(Ax3(k+1,j,i)*Vs(BX3s,k+1,j,i)-Ax3(k,j,i)*Vs(BX3s,k,j,i));  )

      divBmax=FMAX(FABS(D_EXPAND(dB1, +dB2, +dB3))/dV(k,j,i),divBmax);
    },
    Kokkos::Max<real>(divB) // reduction
  );

#ifdef WITH_MPI
  if(idfx::psize>1) {
    MPI_Allreduce(MPI_IN_PLACE, &divB, 1, realMPI, MPI_MAX, MPI_COMM_WORLD);
  }
#endif

  return(divB);
}


/*
real Fluid::CheckDivB(DataBlock &data) {
  real divB=0;
  IdefixArray4D<real> Vs = this->Vs;
  IdefixArray3D<real> Ax1 = data->A[IDIR];
  IdefixArray3D<real> Ax2 = data->A[JDIR];
  IdefixArray3D<real> Ax3 = data->A[KDIR];
  IdefixArray3D<real> dV = data->dV;

  int iref,jref,kref;

  for(int k = data->beg[KDIR] ; k < data->end[KDIR] ; k++) {
    for(int j = data->beg[JDIR] ; j < data->end[JDIR] ; j++) {
      for(int i = data->beg[IDIR] ; i < data->end[IDIR] ; i++) {
        real dB1,dB2,dB3;

        dB1=dB2=dB3=ZERO_F;

        D_EXPAND( dB1=(Ax1(k,j,i+1)*Vs(BX1s,k,j,i+1)-Ax1(k,j,i)*Vs(BX1s,k,j,i)); ,
                  dB2=(Ax2(k,j+1,i)*Vs(BX2s,k,j+1,i)-Ax2(k,j,i)*Vs(BX2s,k,j,i)); ,
                  dB3=(Ax3(k+1,j,i)*Vs(BX3s,k+1,j,i)-Ax3(k,j,i)*Vs(BX3s,k,j,i));  )

        if(FABS(D_EXPAND(dB1, +dB2, +dB3))/dV(k,j,i) > divB) {
          iref=i;
          jref=j;
          kref=k;
          divB=FABS(D_EXPAND(dB1, +dB2, +dB3))/dV(k,j,i);
        }
      }
    }
  }
  idfx::cout << "divB=" << divB << "(i,j,k)=(" << iref << "," << jref << "," << kref
             << ")" << std::endl;
  idfx::cout << "dV=" << dV(kref,jref,iref) << std::endl;
  idfx::cout << " Ax1=" <<Ax1(kref,jref,iref) << " ; " << Ax1(kref,jref,iref+1) << std::endl;
  idfx::cout << " Ax2=" <<Ax2(kref,jref,iref) << " ; " << Ax2(kref,jref+1,iref) << std::endl;
  idfx::cout << " Bx1=" <<Vs(BX1s,kref,jref,iref) << " ; " << Vs(BX1s,kref,jref,iref+1)
             << std::endl;
  idfx::cout << " Bx2=" <<Vs(BX2s,kref,jref,iref) << " ; " << Vs(BX2s,kref,jref+1,iref)
             << std::endl;

  return(divB);
}
*/

#endif // FLUID_CHECKDIVB_HPP_
