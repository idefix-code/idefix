// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>
#include "slice.hpp"
#include "input.hpp"
#include "physics.hpp"
#include "dataBlock.hpp"
#include "fluid.hpp"
#include "vtk.hpp"

Slice::Slice(Input &input, DataBlock & data, int nSlice, SliceType type,
             int direction, real x0, real period) {
  idfx::pushRegion("Slice::Slice");
  std::string prefix = "slice"+std::to_string(nSlice);
  this->slicePeriod = period;
  if(slicePeriod> 0) {
    sliceLast = data.t - slicePeriod;
  }
  // Create the slice.
  this->type = type;
  this->direction = direction;
  // Initialize the subgrid
  this->subgrid = std::make_unique<SubGrid>(data.mygrid, type, direction, x0);
  // Initialize the associated dataBlock
  this->sliceData = std::make_unique<DataBlock>(subgrid.get());
  this->containsX0 = (data.xbeg[direction] <= x0)
                     && (data.xend[direction] >= x0);

  // Initialize the vtk routines
  this->vtk = std::make_unique<Vtk>(input, sliceData.get(),prefix);

  // Allocate array to compute the slice
  this->Vc = IdefixArray4D<real>("Slice_Vc", NVAR,
                                 sliceData->np_tot[KDIR],
                                 sliceData->np_tot[JDIR],
                                 sliceData->np_tot[IDIR]);

  vtk->RegisterVariable(Vc, "RHO", RHO);
  vtk->RegisterVariable(Vc, "VX1", VX1);

  #if MHD == YES
  vtk->RegisterVariable(Vc, "BX1", BX1);
  vtk->RegisterVariable(Vc, "BX2", BX2);
  vtk->RegisterVariable(Vc, "BX3", BX3);
  #endif

  // todo(glesur): add variables for dust and other fluids.

  idfx::popRegion();
}

void Slice::CheckForWrite(DataBlock &data) {
  idfx::pushRegion("Slice:CheckForWrite");

  if(data.t >= sliceLast + slicePeriod) {
    auto Vcin=data.hydro->Vc;
    auto Vcout=this->Vc;
    if(this->type == SliceType::Cut && containsX0) {
      // index of element in current datablock
      const int idx = subgrid->index - data.gbeg[direction]
                                    + data.beg[direction];

      if(direction == IDIR) {
        idefix_for("slice",0,NVAR,0,data.np_tot[KDIR],0,data.np_tot[JDIR],
                  KOKKOS_LAMBDA(int n,int k,int j) {
                    Vcout(n,k,j,0) = Vcin(n,k,j,idx);
                  });
      } else if(direction == JDIR) {
        idefix_for("slice",0,NVAR,0,data.np_tot[KDIR],0,data.np_tot[IDIR],
                  KOKKOS_LAMBDA(int n,int k,int i) {
                    Vcout(n,k,0,i) = Vcin(n,k,idx,i);
                  });
      } else if(direction == KDIR) {
        idefix_for("slice",0,NVAR,0,data.np_tot[JDIR],0,data.np_tot[IDIR],
                  KOKKOS_LAMBDA(int n,int j,int i) {
                    Vcout(n,0,j,i) = Vcin(n,idx,j,i);
                  });
      }
      vtk->Write();
    }
    if(this->type == SliceType::Average) {
      // Perform a point average (NB: this does not perform a volume average!)
      // This is a naive approach which could be optimised using threaded loops
      // However, since this is only for I/O, we don't do any proper optimisation
      idefix_for("Zero",0,NVAR,0,Vcout.extent(1),0,Vcout.extent(2),
        KOKKOS_LAMBDA(int n,int k,int j) {
                    Vcout(n,k,j,0) = 0.0;
                  });
      int beg = data.beg[direction];
      int end = data.end[direction];
      int ntot = data.mygrid->np_int[direction];
      const int dir = direction;
      idefix_for("average",0,NVAR,data.beg[KDIR],data.end[KDIR],
                                  data.beg[JDIR],data.end[JDIR],
                                  data.beg[IDIR],data.end[IDIR],
                  KOKKOS_LAMBDA(int n,int k,int j, int i) {
                    const int it = (dir == IDIR ? 0 : i);
                    const int jt = (dir == JDIR ? 0 : j);
                    const int kt = (dir == KDIR ? 0 : k);

                    Kokkos::atomic_add(&Vcout(n,kt,jt,it) , Vcin(n,k,j,i)/ntot);
                  });
      #ifdef WITH_MPI
        // Create a communicator on which we can do the sum accross processors
        int remainDims[3] = {false, false, false};
        remainDims[direction] = true;
        MPI_Comm avgComm;
        MPI_Cart_sub(subgrid->parentGrid->CartComm, remainDims, &avgComm);
        MPI_Allreduce(MPI_IN_PLACE, Vcout.data(),
                      Vcout.extent(0)*Vcout.extent(1)*Vcout.extent(2)*Vcout.extent(3),
                      realMPI, MPI_SUM, avgComm);
      #endif
      if(containsX0) {
        vtk->Write();
      }
    }

    sliceLast += slicePeriod;
    if((sliceLast+slicePeriod <= data.t) && slicePeriod > 0.0) {
      while(sliceLast <= data.t - slicePeriod) {
        sliceLast += slicePeriod;
      }
    }
  }
  idfx::popRegion();
}
