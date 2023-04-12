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

Slice::Slice(DataBlock & data, int nSlice, SliceType type, int direction, real x0, real period) {
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

  idfx::cout << data.xbeg[direction] <<" "<< data.xend[direction] << std::endl;
  // Initialize the vtk routines
  this->vtk = std::make_unique<Vtk>(sliceData.get(),prefix);

  // Allocate array to compute the slice
  this->Vc = IdefixArray4D<real>("Slice_Vc", NVAR,
                                 sliceData->np_tot[KDIR],
                                 sliceData->np_tot[JDIR],
                                 sliceData->np_tot[IDIR]);

  vtk->RegisterVariable(Vc, "RHO", RHO);
  //vtk->RegisterVariable(Vc, "VX1", VX1);
  idfx::popRegion();
}

void Slice::CheckForWrite(DataBlock &data) {
  idfx::pushRegion("Slice:CheckForWrite");
  // Take the slice (warning with MPI: more work is needed)
  if(data.t >= sliceLast + slicePeriod) {
    idfx::cout << "Go johnny" << std::endl;
    auto Vcin=data.hydro->Vc;
    auto Vcout=this->Vc;
    if(containsX0) {
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

    sliceLast += slicePeriod;
    if((sliceLast+slicePeriod <= data.t) && slicePeriod > 0.0) {
      while(sliceLast <= data.t - slicePeriod) {
        sliceLast += slicePeriod;
      }
    }
  }
  idfx::popRegion();
}
