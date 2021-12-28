// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "dataBlock.hpp"
#include "dataBlockHost.hpp"

// Check if current datablock has nans

int DataBlock::CheckNan()  {
  int nanVs=0;
  int nanVc=0;

  idfx::pushRegion("DataBlock::CheckNan");
  IdefixArray4D<real> Vc=this->hydro.Vc;

  idefix_reduce("checkNanVc",
    0, NVAR,
    beg[KDIR], end[KDIR],
    beg[JDIR], end[JDIR],
    beg[IDIR], end[IDIR],
    KOKKOS_LAMBDA (int n, int k, int j, int i, int &nnan) {
      if(std::isnan(Vc(n,k,j,i))) nnan++;
    }, Kokkos::Sum<int>(nanVc) // reduction variable
  );

  #if MHD == YES
    IdefixArray4D<real> Vs=this->hydro.Vs;
    idefix_reduce("checkNanVs",
      0, DIMENSIONS,
      beg[KDIR], end[KDIR]+KOFFSET,
      beg[JDIR], end[JDIR]+JOFFSET,
      beg[IDIR], end[IDIR]+IOFFSET,
      KOKKOS_LAMBDA (int n, int k, int j, int i, int &nnan) {
        if(std::isnan(Vs(n,k,j,i))) nnan++;
      }, Kokkos::Sum<int>(nanVs) // reduction variable
    );
  #endif

  int nanTot = nanVc+nanVs;
  #ifdef WITH_MPI
    MPI_Allreduce(MPI_IN_PLACE, &nanTot,1,MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  #endif
  if(nanTot>0) {
    idfx::cout << "DataBlock: Nans were found in the current calculation" << std::endl;
    if(nanVc+nanVs>0) {
      // Each datablock shows its findings to stdout
      std::cout << "DataBlock: rank " << idfx::prank << " found " << nanVc << " Nans in Vc"
      #if MHD == YES
        << " and " << nanVs << " Nans in Vs"
      #endif
        << " in the current datablock. Details will be in corresponding process log file."
        << std::endl;

      // We need to make copies to find exactly where the thing is wrong
      DataBlockHost dataHost(*this);
      dataHost.SyncFromDevice();

      int nerrormax=10;

      for(int k = beg[KDIR] ; k < end[KDIR] ; k++) {
        for(int j = beg[JDIR] ; j < end[JDIR] ; j++) {
          for(int i = beg[IDIR] ; i < end[IDIR] ; i++) {
            for(int n = 0 ; n < NVAR ; n ++) {
              if(std::isnan(dataHost.Vc(n,k,j,i)) && nerrormax>0) {
                nerrormax--;
                idfx::cout << "rank " << idfx::prank << ": Nan found  in variable "
                  << this->hydro.VcName[n] << std::endl;

                idfx::cout << "      global (i,j,k) = (" << i-beg[IDIR]+gbeg[IDIR]-nghost[IDIR]
                  << ", " << j-beg[JDIR]+gbeg[JDIR]-nghost[JDIR] << ", "
                  << k-beg[KDIR]+gbeg[KDIR]-nghost[KDIR] << ")" << std::endl;

                idfx::cout << "      global (x,y,z) = (" << dataHost.x[IDIR](i) << ", "
                  <<  dataHost.x[JDIR](j) << ", " << dataHost.x[KDIR](k) << ")" << std::endl;
              }
            }
          }
        }
      }

  #if MHD == YES
      for(int k = beg[KDIR] ; k < end[KDIR]+KOFFSET ; k++) {
        for(int j = beg[JDIR] ; j < end[JDIR]+JOFFSET ; j++) {
          for(int i = beg[IDIR] ; i < end[IDIR]+IOFFSET ; i++) {
            for(int n = 0 ; n < DIMENSIONS ; n ++) {
              if(std::isnan(dataHost.Vs(n,k,j,i)) && nerrormax>0) {
                nerrormax--;
                idfx::cout << "rank " << idfx::prank << ": Nan found  in variable "
                  << this->hydro.VsName[n] << std::endl;
                idfx::cout << "      global (i,j,k) = (" << i-beg[IDIR]+gbeg[IDIR]-nghost[IDIR]
                  << ", " << j-beg[JDIR]+gbeg[JDIR]-nghost[JDIR] << ", "
                  << k-beg[KDIR]+gbeg[KDIR]-nghost[KDIR] << ")" << std::endl;

                idfx::cout << "      global (x,y,z) = (" << dataHost.x[IDIR](i) << ", "
                  <<  dataHost.x[JDIR](j) << ", " << dataHost.x[KDIR](k) << ")" << std::endl;
              }
            }
          }
        }
      }
  #endif
      if(nerrormax<=0) {
        idfx::cout << "... " << std::endl << "*** More Nans have been found in current dataBlock. "
          << "Only showing the first 10." << std::endl;
      }
    }
  }
  idfx::popRegion();
  return(nanTot);
}
