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

  // Allocate array to compute the slice
  this->Vc = IdefixArray4D<real>("Slice_Vc", NVAR,
                                 sliceData->np_tot[KDIR],
                                 sliceData->np_tot[JDIR],
                                 sliceData->np_tot[IDIR]);

  vtk->RegisterVariable(Vc, "RHO", RHO);
  EXPAND( vtk->RegisterVariable(Vc, "VX1", VX1);   ,
          vtk->RegisterVariable(Vc, "VX2", VX2);   ,
          vtk->RegisterVariable(Vc, "VX3", VX3); )
  #if MHD == YES
  EXPAND( vtk->RegisterVariable(Vc, "BX1", BX1);   ,
          vtk->RegisterVariable(Vc, "BX2", BX2);   ,
          vtk->RegisterVariable(Vc, "BX3", BX3); )
  #endif

  // todo(glesur): add variables for dust and other fluids.

  idfx::popRegion();
}

void Slice::EnrollUserDefVariables(std::map<std::string,IdefixHostArray3D<real>> userDefVar) {
  this->userDefVariableFull = userDefVar;
  this->haveUserDefinedVariables = true;
  for(auto const &[name, array] : userDefVariableFull) {
    userDefVariableSliced.emplace(name, IdefixHostArray3D<real>(std::string("Slice_"+name),
                                                                sliceData->np_tot[KDIR],
                                                                sliceData->np_tot[JDIR],
                                                                sliceData->np_tot[IDIR]));
    vtk->RegisterVariable(userDefVariableSliced[name],name);
  }
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
      userDefVariablesFunc(data, userDefVariableFull);
      idfx::popRegion();
    }
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
        // Take a slice of the user-defined variables
        if(haveUserDefinedVariables) {
          for(auto const &[name, array] : userDefVariableFull) {
            for(int k = 0 ; k < data.np_tot[KDIR] ; k++) {
              for(int j = 0 ; j < data.np_tot[JDIR] ; j++) {
                userDefVariableSliced[name](k,j,0) = userDefVariableFull[name](k,j,idx);
              }
            }
          }
        }
      } else if(direction == JDIR) {
        idefix_for("slice",0,NVAR,0,data.np_tot[KDIR],0,data.np_tot[IDIR],
                  KOKKOS_LAMBDA(int n,int k,int i) {
                    Vcout(n,k,0,i) = Vcin(n,k,idx,i);
                  });
        // Take a slice of the user-defined variables
        if(haveUserDefinedVariables) {
          for(auto const &[name, array] : userDefVariableFull) {
            for(int k = 0 ; k < data.np_tot[KDIR] ; k++) {
              for(int i = 0 ; i < data.np_tot[IDIR] ; i++) {
                userDefVariableSliced[name](k,0,i) = userDefVariableFull[name](k,idx,i);
              }
            }
          }
        }
      } else if(direction == KDIR) {
        idefix_for("slice",0,NVAR,0,data.np_tot[JDIR],0,data.np_tot[IDIR],
                  KOKKOS_LAMBDA(int n,int j,int i) {
                    Vcout(n,0,j,i) = Vcin(n,idx,j,i);
                  });
        // Take a slice of the user-defined variables
        if(haveUserDefinedVariables) {
          for(auto const &[name, array] : userDefVariableFull) {
            for(int j = 0 ; j < data.np_tot[JDIR] ; j++) {
              for(int i = 0 ; i < data.np_tot[IDIR] ; i++) {
                userDefVariableSliced[name](0,j,i) = userDefVariableFull[name](idx,j,i);
              }
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
      idefix_for("Zero",0,NVAR,0,Vcout.extent(1),0,Vcout.extent(2),0,Vcout.extent(3),
        KOKKOS_LAMBDA(int n,int k,int j, int i) {
                    Vcout(n,k,j,i) = 0.0;
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
        Kokkos::fence();
        MPI_Allreduce(MPI_IN_PLACE, Vcout.data(),
                      Vcout.extent(0)*Vcout.extent(1)*Vcout.extent(2)*Vcout.extent(3),
                      realMPI, MPI_SUM, avgComm);
      #endif
      // Average of user-defined variables
      if(haveUserDefinedVariables) {
        for(auto const &[name, array] : userDefVariableSliced) {
          for(int k = 0 ; k < array.extent(0) ; k++) {
            for(int j = 0 ; j < array.extent(1) ; j++) {
              for(int i = 0 ; i < array.extent(2) ; i++) {
                array(k,j,i) = 0.0;
          }}}
          for(int k = data.beg[KDIR] ; k < data.end[KDIR] ; k++) {
            for(int j = data.beg[JDIR] ; j < data.end[JDIR] ; j++) {
              for(int i = data.beg[IDIR] ; i < data.end[IDIR] ; i++) {
                const int it = (dir == IDIR ? 0 : i);
                const int jt = (dir == JDIR ? 0 : j);
                const int kt = (dir == KDIR ? 0 : k);
                array(kt,jt,it) += userDefVariableFull[name](k,j,i)/ntot;
          }}}
          #ifdef WITH_MPI
          MPI_Allreduce(MPI_IN_PLACE, array.data(),
                    array.extent(0)*array.extent(1)*array.extent(2),
                    realMPI, MPI_SUM, avgComm);
          #endif
        }
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
