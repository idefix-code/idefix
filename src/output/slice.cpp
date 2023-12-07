// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>
#include <map>
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
  // Register the last output in dumps so that we restart from the right slice
  data.dump->RegisterVariable(&sliceLast, std::string("slcLast-")+std::to_string(nSlice));
  // Create the slice.
  this->type = type;
  this->direction = direction;
  // Initialize the subgrid
  this->subgrid = std::make_unique<SubGrid>(data.mygrid, type, direction, x0);
  // Initialize the associated dataBlock
  this->sliceData = std::make_unique<DataBlock>(subgrid.get());
  this->containsX0 = (data.xbeg[direction] <= x0)
                     && (data.xend[direction] > x0);

  #ifdef WITH_MPI
    if(type==SliceType::Average) {
      // Create a communicator on which we can do the sum accross processors
      int remainDims[3] = {false, false, false};
      remainDims[direction] = true;
      MPI_Cart_sub(subgrid->parentGrid->CartComm, remainDims, &avgComm);
    }
  #endif


  // Initialize the vtk routines
  this->vtk = std::make_unique<Vtk>(input, sliceData.get(),prefix);
  // Make sure the vtk file gets store in the parent's datablock dump
  data.dump->RegisterVariable(&vtk->vtkFileNumber,
                              std::string("slcNumber-")+std::to_string(nSlice));


  // Allocate array to compute the slice of each variable registered for VTK output
  // in the parent dataBlock
  for(auto const& [name, scalar] : data.vtk->vtkScalarMap) {
    IdefixHostArray3D<real> arr("Slice_Vc",
                                 sliceData->np_tot[KDIR],
                                 sliceData->np_tot[JDIR],
                                 sliceData->np_tot[IDIR]);
    this->variableMap.emplace(name, arr);
    vtk->RegisterVariable(arr,name);
  }
  // todo(glesur): add variables for dust and other fluids.

  idfx::popRegion();
}

void Slice::EnrollUserDefVariables(std::map<std::string,IdefixHostArray3D<real>> userDefVar) {
  this->userDefVariableMap = userDefVar;
  this->haveUserDefinedVariables = true;
}

void Slice::EnrollUserDefFunc(UserDefVariablesFunc myFunc) {
  userDefVariablesFunc = myFunc;
}

void Slice::CheckForWrite(DataBlock &data) {
  idfx::pushRegion("Slice:CheckForWrite");

  if(data.t >= sliceLast + slicePeriod) {
    // sync time
    sliceData->t = data.t;
    if(haveUserDefinedVariables) {
      // Call user-def function to fill the userdefined variable arrays
      idfx::pushRegion("UserDef::User-defined variables function");
      userDefVariablesFunc(data, userDefVariableMap);
      idfx::popRegion();
    }

    if(this->type == SliceType::Cut && containsX0) {
      // index of element in current datablock
      const int idx = subgrid->index - data.gbeg[direction]
                                    + data.beg[direction];

      if(direction == IDIR) {
        for(auto const &it : variableMap) {
          auto name = it.first;
          auto &scalar= data.vtk->vtkScalarMap.find(name)->second;
          auto arrIn = scalar.GetHostField();
          auto arrOut = variableMap[name];
          for(int k = 0 ; k < data.np_tot[KDIR] ; k++) {
            for(int j = 0 ; j < data.np_tot[JDIR] ; j++) {
              arrOut(k,j,0) = arrIn(k,j,idx);
            }
          }
        }
      } else if(direction == JDIR) {
        for(auto const &it : variableMap) {
          auto name = it.first;
          auto &scalar= data.vtk->vtkScalarMap.find(name)->second;
          auto arrIn = scalar.GetHostField();
          auto arrOut = variableMap[name];
          for(int k = 0 ; k < data.np_tot[KDIR] ; k++) {
            for(int i = 0 ; i < data.np_tot[IDIR] ; i++) {
              arrOut(k,0,i) = arrIn(k,idx,i);
            }
          }
        }
      } else if(direction == KDIR) {
        for(auto const &it : variableMap) {
          auto name = it.first;
          auto &scalar= data.vtk->vtkScalarMap.find(name)->second;
          auto arrIn = scalar.GetHostField();
          auto arrOut = variableMap[name];
          for(int j = 0 ; j < data.np_tot[JDIR] ; j++) {
            for(int i = 0 ; i < data.np_tot[IDIR] ; i++) {
              arrOut(0,j,i) = arrIn(idx,j,i);
            }
          }
        }
      }
      vtk->Write();
    }
    if(this->type == SliceType::Average) {
      // Perform a point average (NB: this does not perform a volume average!)
      // This is a naive approach which could be optimised using threaded loops
      // However, since this is only for I/O, we don't do any proper optimisation
      int beg = data.beg[direction];
      int end = data.end[direction];
      int ntot = data.mygrid->np_int[direction];
      const int dir = direction;
      // Averaged of registered varaibles
      for(auto const &it : variableMap) {
        auto name = it.first;
        auto &scalar= data.vtk->vtkScalarMap.find(name)->second;
        auto arrIn = scalar.GetHostField();
        auto arrOut = variableMap[name];
        for(int k = 0 ; k < arrOut.extent(0) ; k++) {
          for(int j = 0 ; j < arrOut.extent(1) ; j++) {
            for(int i = 0 ; i < arrOut.extent(2) ; i++) {
              arrOut(k,j,i) = 0.0;
        }}}
        for(int k = data.beg[KDIR] ; k < data.end[KDIR] ; k++) {
          for(int j = data.beg[JDIR] ; j < data.end[JDIR] ; j++) {
            for(int i = data.beg[IDIR] ; i < data.end[IDIR] ; i++) {
              const int it = (dir == IDIR ? 0 : i);
              const int jt = (dir == JDIR ? 0 : j);
              const int kt = (dir == KDIR ? 0 : k);
              arrOut(kt,jt,it) += arrIn(k,j,i)/ntot;
        }}}
        #ifdef WITH_MPI
          Kokkos::fence();
          MPI_Allreduce(MPI_IN_PLACE, arrOut.data(),
                        arrOut.extent(0)*arrOut.extent(1)*arrOut.extent(2)*arrOut.extent(3),
                        realMPI, MPI_SUM, avgComm);
        #endif
      }
      if(containsX0) {
        vtk->Write();
      }
    }
    if(!containsX0) {
      vtk->vtkFileNumber++; // increment file number so that each process stay in sync
    }

    sliceLast += slicePeriod;
    if((sliceLast+slicePeriod <= data.t) && slicePeriod > 0.0) {
      while(sliceLast <= data.t - slicePeriod) {
        sliceLast += slicePeriod;
      }
    }
    #ifdef WITH_MPI
      MPI_Barrier(MPI_COMM_WORLD);
    #endif
  }
  idfx::popRegion();
}
