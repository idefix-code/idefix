// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_CHECKNAN_HPP_
#define FLUID_CHECKNAN_HPP_
#include <iostream>
#include "dataBlock.hpp"
#include "dataBlockHost.hpp"
#include "fluid.hpp"

// Check if current datablock has nans

template<typename Phys>
int Fluid<Phys>::CheckNan()  {
  int nanVs=0;
  int nanVc=0;

  idfx::pushRegion("Fluid::CheckNan");
  IdefixArray4D<real> Vc=this->Vc;

  idefix_reduce("checkNanVc",
    0, Phys::nvar,
    data->beg[KDIR], data->end[KDIR],
    data->beg[JDIR], data->end[JDIR],
    data->beg[IDIR], data->end[IDIR],
    KOKKOS_LAMBDA (int n, int k, int j, int i, int &nnan) {
      if(std::isnan(Vc(n,k,j,i))) nnan++;
    }, Kokkos::Sum<int>(nanVc) // reduction variable
  );

  if constexpr(Phys::mhd) {
    IdefixArray4D<real> Vs=this->Vs;
    idefix_reduce("checkNanVs",
      0, DIMENSIONS,
      data->beg[KDIR], data->end[KDIR]+KOFFSET,
      data->beg[JDIR], data->end[JDIR]+JOFFSET,
      data->beg[IDIR], data->end[IDIR]+IOFFSET,
      KOKKOS_LAMBDA (int n, int k, int j, int i, int &nnan) {
        if(std::isnan(Vs(n,k,j,i))) nnan++;
      }, Kokkos::Sum<int>(nanVs) // reduction variable
    );
  }


  int nanTot = nanVc+nanVs;
  #ifdef WITH_MPI
    MPI_Allreduce(MPI_IN_PLACE, &nanTot,1,MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  #endif
  if(nanTot>0) {
    idfx::cout << "Fluid<" << prefix << ">: Nans were found in the current calculation"
               << std::endl;
    if(nanVc+nanVs>0) {
      // Each datablock shows its findings to stdout
      std::cout << "Fluid<" << prefix << ">: rank " << idfx::prank
                << " found " << nanVc << " Nans in Vc"
      #if MHD == YES
        << " and " << nanVs << " Nans in Vs"
      #endif
        << " in the current datablock. Details will be in corresponding process log file."
        << std::endl;

      // We need to make copies to find exactly where the thing is wrong

      DataBlockHost dataHost(*data);

      IdefixHostArray4D<real> VcHost = Kokkos::create_mirror_view(this->Vc);
      Kokkos::deep_copy(VcHost,Vc);

      int nerrormax=10;

      for(int k = data->beg[KDIR] ; k < data->end[KDIR] ; k++) {
        for(int j = data->beg[JDIR] ; j < data->end[JDIR] ; j++) {
          for(int i = data->beg[IDIR] ; i < data->end[IDIR] ; i++) {
            for(int n = 0 ; n < Phys::nvar ; n ++) {
              if(std::isnan(VcHost(n,k,j,i)) && nerrormax>0) {
                nerrormax--;
                idfx::cout << "rank " << idfx::prank << ": Nan found  in variable "
                  << this->VcName[n] << std::endl;

                idfx::cout << "      global (i,j,k) = ("
                  << i-data->beg[IDIR]+data->gbeg[IDIR]-data->nghost[IDIR]
                  << ", " << j-data->beg[JDIR]+data->gbeg[JDIR]-data->nghost[JDIR] << ", "
                  << k-data->beg[KDIR]+data->gbeg[KDIR]-data->nghost[KDIR] << ")" << std::endl;

                idfx::cout << "      global (x,y,z) = (" << dataHost.x[IDIR](i) << ", "
                  <<  dataHost.x[JDIR](j) << ", " << dataHost.x[KDIR](k) << ")" << std::endl;
              }
            }
          }
        }
      }

      if constexpr(Phys::mhd) {
        IdefixHostArray4D<real> VsHost = Kokkos::create_mirror_view(this->Vs);
        Kokkos::deep_copy(VsHost,Vs);
        for(int k = data->beg[KDIR] ; k < data->end[KDIR]+KOFFSET ; k++) {
          for(int j = data->beg[JDIR] ; j < data->end[JDIR]+JOFFSET ; j++) {
            for(int i = data->beg[IDIR] ; i < data->end[IDIR]+IOFFSET ; i++) {
              for(int n = 0 ; n < DIMENSIONS ; n ++) {
                if(std::isnan(VsHost(n,k,j,i)) && nerrormax>0) {
                  nerrormax--;
                  idfx::cout << "rank " << idfx::prank << ": Nan found  in variable "
                    << this->VsName[n] << std::endl;
                  idfx::cout << "      global (i,j,k) = ("
                             << i-data->beg[IDIR]+data->gbeg[IDIR]-data->nghost[IDIR]
                    << ", " << j-data->beg[JDIR]+data->gbeg[JDIR]-data->nghost[JDIR] << ", "
                    << k-data->beg[KDIR]+data->gbeg[KDIR]-data->nghost[KDIR] << ")" << std::endl;

                  idfx::cout << "      global (x,y,z) = (" << dataHost.x[IDIR](i) << ", "
                    <<  dataHost.x[JDIR](j) << ", " << dataHost.x[KDIR](k) << ")" << std::endl;
                }
              }
            }
          }
        }
      }
      if(nerrormax<=0) {
        idfx::cout << "... " << std::endl << "*** More Nans have been found in current fluid. "
          << "Only showing the first 10." << std::endl;
      }
    }
  }
  idfx::popRegion();
  return(nanTot);
}

#endif // FLUID_CHECKNAN_HPP_
