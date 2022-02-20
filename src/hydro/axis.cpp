// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <vector>
#include "axis.hpp"
#include "hydro.hpp"
#include "dataBlock.hpp"


void Axis::Init(Grid &grid, Hydro *h) {
  this->hydro = h;
  this->data = this->hydro->data;
  this->emf = & this->hydro->emf;

  #if GEOMETRY != SPHERICAL
    IDEFIX_ERROR("Axis boundary conditions are only designed to handle spherical geometry");
  #endif
  idfx::cout << "Axis: Axis regularisation enabled ";

  if(fabs((grid.xend[KDIR] - grid.xbeg[KDIR] -2.0*M_PI)) < 1e-10) {
    this->isTwoPi = true;
    idfx::cout << "with full (2pi) azimuthal extension" << std::endl;;
    #ifdef WITH_MPI
      // Check that there is a domain decomposition in phi
      if(data->mygrid->nproc[KDIR]>1) {
        if(data->mygrid->nproc[KDIR]%2==1) {
          IDEFIX_ERROR("The numbre of processes in the phi direction should"
                        " be even for axis decomposition");
        }
        needMPIExchange = true;
      }
    #endif
  } else {
    this->isTwoPi = false;
    idfx::cout << "with partial (<2pi) azimuthal extension" << std::endl;
  }

  // Check where the axis is lying.
  if(hydro->data->lbound[JDIR] == axis) axisLeft = true;
  if(hydro->data->rbound[JDIR] == axis) axisRight = true;

  // Init the symmetry array (used to flip the signs of arrays accross the axis)
  symmetryVc = IdefixArray1D<int>("Axis:SymmetryVc",NVAR);
  IdefixArray1D<int>::HostMirror symmetryVcHost = Kokkos::create_mirror_view(symmetryVc);
  // Init the array
  for (int nv = 0; nv < NVAR; nv++) {
    symmetryVcHost(nv) = 1;
    if (nv == VX2)
      symmetryVcHost(nv) = -1;
    if (nv == VX3)
      symmetryVcHost(nv) = -1;
    if (nv == BX2)
      symmetryVcHost(nv) = -1;
    if (nv == BX3)
      symmetryVcHost(nv) = -1;
  }
  Kokkos::deep_copy(symmetryVc, symmetryVcHost);

#if MHD == YES
  symmetryVs = IdefixArray1D<int>("Axis:SymmetryVs",DIMENSIONS);
  IdefixArray1D<int>::HostMirror symmetryVsHost = Kokkos::create_mirror_view(symmetryVs);
  // Init the array
  for(int nv = 0; nv < DIMENSIONS; nv++) {
    symmetryVsHost(nv) = 1;
    if (nv == BX2s)
      symmetryVsHost(nv) = -1;
    if (nv == BX3s)
      symmetryVsHost(nv) = -1;
  }
  Kokkos::deep_copy(symmetryVs, symmetryVsHost);
#endif

  #if MHD == YES
  this->Ex1Avg = IdefixArray1D<real>("Axis:Ex1Avg",hydro->data->np_tot[IDIR]);
  #endif

  #ifdef WITH_MPI
    if(needMPIExchange) {
      // Make MPI exchange datatypes
      MakeMPIDataypes(JDIR);
    }
  #endif
}

void Axis::SymmetrizeEx1Side(int jref) {
#if DIMENSIONS == 3
  IdefixArray3D<real> Ex1 = emf->ex;
  IdefixAtomicArray1D<real> Ex1Avg = this->Ex1Avg;

  if(isTwoPi) {
    idefix_for("Ex1_ini",0,data->np_tot[IDIR],
        KOKKOS_LAMBDA(int i) {
          Ex1Avg(i) = ZERO_F;
        });

    idefix_for("Ex1_Symmetrize",data->beg[KDIR],data->end[KDIR],0,data->np_tot[IDIR],
      KOKKOS_LAMBDA(int k,int i) {
        Ex1Avg(i) += Ex1(k,jref,i);
      });
    if(needMPIExchange) {
      #ifdef WITH_MPI
        // sum along all of the processes on the same r
        MPI_Allreduce(MPI_IN_PLACE, Ex1Avg.data(), data->np_tot[IDIR], realMPI,
                      MPI_SUM, data->mygrid->AxisComm);
      #endif
    }

    int ncells=data->mygrid->np_int[KDIR];

    idefix_for("Ex1_Store",0,data->np_tot[KDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA(int k,int i) {
      Ex1(k,jref,i) = Ex1Avg(i)/((real) ncells);
    });
  } else {
    // if we're not doing full two pi, the flow is symmetric with respect to the axis, and the axis
    // EMF is simply zero
    idefix_for("Ex1_Store",0,data->np_tot[KDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA(int k,int i) {
      Ex1(k,jref,i) = ZERO_F;
    });
  }
#endif
}

// Average the Emf component along the axis
void Axis::SymmetrizeEx1() {
  idfx::pushRegion("Axis::SymmetrizeEx1");

  if(this->axisLeft) {
    int jref = hydro->data->beg[JDIR];
    this->SymmetrizeEx1Side(jref);
  }
  if(this->axisRight) {
    int jref = hydro->data->end[JDIR];
    this->SymmetrizeEx1Side(jref);
  }

  idfx::popRegion();
}

// enforce the boundary conditions on the ghost zone accross the axis
void Axis::EnforceAxisBoundary(int side) {
  idfx::pushRegion("Axis::EnforceAxisBoundary");
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray1D<int> sVc = this->symmetryVc;

  int ibeg = 0;
  int iend = data->np_tot[IDIR];
  int jref, jbeg,jend;
  int offset;
  if(side == left) {
    jref = data->beg[JDIR];
    jbeg = 0;
    jend = data->beg[JDIR];
    offset = -1;
  }
  if(side == right) {
    jref = data->end[JDIR]-1;
    jbeg = data->end[JDIR];
    jend = data->np_tot[JDIR];
    offset = 1;
  }

  int kbeg = 0;
  int kend = data->np_tot[KDIR];
  int np_int_k = data->mygrid->np_int[KDIR];
  int nghost_k = data->mygrid->nghost[KDIR];

  if(isTwoPi) {
    if(needMPIExchange) {
      ExchangeMPI(side);
      idefix_for("BoundaryAxis",0,NVAR,kbeg,kend,ibeg,iend,
        KOKKOS_LAMBDA (int n, int k, int i) {
          real scrch[4];    // scratch array (max 4 elements)
          // Exchange axis copy the required array, but we must invert the JDIR axis
          for(int j = jbeg ; j < jend ; j++) {
            scrch[j-jbeg] =  Vc(n,k,j,i);
          }
          for(int j = jbeg ; j < jend ; j++) {
            Vc(n,k,j,i) = sVc(n) * scrch[(jend-1) - j];
          }
        });
    } else { // no MPI exchange
      idefix_for("BoundaryAxis",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
              KOKKOS_LAMBDA (int n, int k, int j, int i) {
                int kcomp = nghost_k + (( k - nghost_k + np_int_k/2) % np_int_k);

                Vc(n,k,j,i) = sVc(n)*Vc(n, kcomp, 2*jref-j+offset,i);
              });
    }// MPI Exchange
  } else {  // not 2pi
    idefix_for("BoundaryAxis",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
            KOKKOS_LAMBDA (int n, int k, int j, int i) {
              // kcomp = k by construction since we're doing a fraction of twopi

              Vc(n,k,j,i) = sVc(n)*Vc(n, k, 2*jref-j+offset,i);
            });
  }

  #if MHD == YES
    IdefixArray4D<real> Vs = hydro->Vs;
    IdefixArray1D<int> sVs = this->symmetryVs;

    for(int component=0; component<DIMENSIONS; component++) {
      int ieb,jeb,keb;
      if(component == IDIR)
        ieb=iend+1;
      else
        ieb=iend;
      if(component == JDIR)
        jeb=jend+1;
      else
        jeb=jend;
      if(component == KDIR)
        keb=kend+1;
      else
        keb=kend;
      if(component != JDIR) { // skip normal component
        if(isTwoPi) {
          if(needMPIExchange) {
            // MPI exchange already performed during Vc
            idefix_for("BoundaryAxisVs",kbeg,keb,ibeg,ieb,
              KOKKOS_LAMBDA (int k, int i) {
                real scrch[4];    // scratch array (max 4 elements)
                // Exchange axis copy the required array, but we must invert the JDIR axis
                for(int j = jbeg ; j < jend ; j++) {
                  scrch[j-jbeg] =  Vs(component,k,j,i);
                }
                for(int j = jbeg ; j < jend ; j++) {
                  Vs(component,k,j,i) = sVs(component) * scrch[(jend-1) - j];
                }
              });
          } else { //no mpi exchange
            idefix_for("BoundaryAxisVs",kbeg,keb,jbeg,jeb,ibeg,ieb,
              KOKKOS_LAMBDA (int k, int j, int i) {
                int kcomp = nghost_k + (( k - nghost_k + np_int_k/2) % np_int_k);
                Vs(component,k,j,i) = sVs(component)*Vs(component,kcomp, 2*jref-j+offset,i);
              }
            );
          } // mpi exchange
        } else { // not 2pi
          idefix_for("BoundaryAxisVs",kbeg,keb,jbeg,jeb,ibeg,ieb,
            KOKKOS_LAMBDA (int k, int j, int i) {
              Vs(component,k,j,i) = sVs(component)*Vs(component,k, 2*jref-j+offset,i);
            }
          );
        }
      }
    }
  #endif


  idfx::popRegion();
}

// Reconstruct Bx2s taking care of the sides where an axis is lying
void Axis::ReconstructBx2s() {
  idfx::pushRegion("Axis::ReconstructBx2s");
  IdefixArray4D<real> Vs = hydro->Vs;
  IdefixArray3D<real> Ax1=data->A[IDIR];
  IdefixArray3D<real> Ax2=data->A[JDIR];
  IdefixArray3D<real> Ax3=data->A[KDIR];
  int nstart = data->beg[JDIR]-1;
  int nend = data->end[JDIR];
  int ntot = data->np_tot[JDIR];

  int signLeft = 1;
  int signRight = 1;
  if(axisLeft) signLeft = -1;
  if(axisRight) signRight = -1;

#if DIMENSIONS >= 2
  // This loop is a copy of ReconstructNormalField, with the proper sign when we cross the axis
  idefix_for("Axis::ReconstructBX2s",0,data->np_tot[KDIR],0,data->np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int i) {
          for(int j = nstart ; j>=0 ; j-- ) {
            Vs(BX2s,k,j,i) = 1.0 / Ax2(k,j,i) * ( Ax2(k,j+1,i)*Vs(BX2s,k,j+1,i)
                        +(D_EXPAND( Ax1(k,j,i+1) * Vs(BX1s,k,j,i+1) - Ax1(k,j,i) * Vs(BX1s,k,j,i)  ,
                                                                                                  ,
                              + signLeft*( Ax3(k+1,j,i) * Vs(BX3s,k+1,j,i)
                                                              - Ax3(k,j,i) * Vs(BX3s,k,j,i) ))));
          }
          for(int j = nend ; j<ntot ; j++ ) {
            Vs(BX2s,k,j+1,i) = 1.0 / Ax2(k,j+1,i) * ( Ax2(k,j,i)*Vs(BX2s,k,j,i)
                        -(D_EXPAND( Ax1(k,j,i+1) * Vs(BX1s,k,j,i+1) - Ax1(k,j,i) * Vs(BX1s,k,j,i)  ,
                                                                                                  ,
                              + signRight*( Ax3(k+1,j,i) * Vs(BX3s,k+1,j,i)
                                                               - Ax3(k,j,i) * Vs(BX3s,k,j,i) ))));
          }
        }
      );

  // Set BX2s on the axis to the average of the two agacent cells
  // This is required since Bx2s on the axis is not evolved since
  // there is no circulation around it
    bool left = axisLeft;
    bool right = axisRight;

    int jright = data->end[JDIR];
    int jleft = data->beg[JDIR];
    if(isTwoPi) {
      idefix_for("Axis:BoundaryAvg",0,data->np_tot[KDIR],0,data->np_tot[IDIR],
            KOKKOS_LAMBDA (int k, int i) {
              if(left) {
                Vs(BX2s,k,jleft,i) = HALF_F*(Vs(BX2s,k,jleft-1,i)+Vs(BX2s,k,jleft+1,i));
              }
              if(right) {
                Vs(BX2s,k,jright,i) = HALF_F*(Vs(BX2s,k,jright-1,i)+Vs(BX2s,k,jright+1,i));
              }
            }
          );
    } else {
      idefix_for("Axis:BoundaryAvg",0,data->np_tot[KDIR],0,data->np_tot[IDIR],
            KOKKOS_LAMBDA (int k, int i) {
              if(left) {
                Vs(BX2s,k,jleft,i) = ZERO_F;
              }
              if(right) {
                Vs(BX2s,k,jright,i) = ZERO_F;
              }
            }
          );
    }
#endif


  idfx::popRegion();
}

void Axis::ExchangeMPI(int side) {
  idfx::pushRegion("Axis::ExchangeMPI");

#ifdef WITH_MPI
  idfx::mpiCallsTimer -= MPI_Wtime();
  int procSend, procRecv;

  std::vector<MPI_Request> request;


  // We receive from procRecv, and we send to procSend, send to the right. Shift by half the domain
  MPI_SAFE_CALL(MPI_Cart_shift(data->mygrid->AxisComm,0,data->mygrid->nproc[KDIR]/2,
                               &procRecv,&procSend ));

  request.emplace_back();
  MPI_SAFE_CALL(MPI_Irecv(hydro->Vc.data(), 1, typeVcRecv[side], procRecv,
                 9000+10*side, data->mygrid->AxisComm, &request.back()));

  request.emplace_back();
  MPI_SAFE_CALL(MPI_Isend(hydro->Vc.data(), 1, typeVcSend[side], procSend,
                9000+10*side, data->mygrid->AxisComm, &request.back()));

  #if MHD==YES
    request.emplace_back();
    MPI_SAFE_CALL(MPI_Irecv(hydro->Vs.data(), 1, typeVsRecv[side], procRecv,
                  9000+10*side+5, data->mygrid->AxisComm, &request.back()));

    request.emplace_back();
    MPI_SAFE_CALL(MPI_Isend(hydro->Vs.data(), 1, typeVsSend[side], procSend,
                  9000+10*side+5, data->mygrid->AxisComm, &request.back()));
  #endif


  std::vector<MPI_Status> status(request.size());
  MPI_Waitall(request.size(), request.data(), status.data());


  idfx::mpiCallsTimer += MPI_Wtime();
#endif
  idfx::popRegion();
}

void Axis::MakeMPIDataypes(int dir) {
  idfx::pushRegion("Axis::MakeMPIDataypes");
#ifdef WITH_MPI
  // Init the datatype in each direction
  int size[3];
  int startSend[3];
  int startRecv[3];
  int subsize[3];

  ////////////////////////////////////////////////////////////////////////////
  // Init variable mappers
  // The variable mapper list all of the variable which are exchanged in MPI boundary calls
  // This is required since we skip some of the variables in Vc to limit the amount of data
  // being exchanged
  #if MHD == YES
  int mapNVars = NVAR - DIMENSIONS; // We will not send magnetic field components which are in Vs
  #else
  int mapNVars = NVAR;
  #endif

  std::vector<int> mapVars;
  // Init the list of variables we will exchange in MPI routines
  int ntarget = 0;
  for(int n = 0 ; n < mapNVars ; n++) {
    mapVars.push_back(ntarget);
    ntarget++;
    #if MHD == YES
      // Skip centered field components if they are also defined in Vs
      #if DIMENSIONS >= 1
        if(ntarget==BX1) ntarget++;
      #endif
      #if DIMENSIONS >= 2
        if(ntarget==BX2) ntarget++;
      #endif
      #if DIMENSIONS == 3
        if(ntarget==BX3) ntarget++;
      #endif
    #endif
  }

  ///////////////////////////////////////////////////////////////////////

  // Create two communicators (one for each face)
  typeVcSend = std::vector<MPI_Datatype>(2);
  typeVcRecv = std::vector<MPI_Datatype>(2);
  // the direction of exchange is assumed to be JDIR
  for(int face = 0 ; face <= 1 ; face++) {
    // We first define the sub-array for one single variable
    for(int i = 0 ; i < 3 ; i++) {
      int ic = 2-i; // because X1 is actually the last index in our arrays
      size[ic] = data->np_tot[i];
      if(i<dir) {
        startRecv[ic] = startSend[ic] = 0;
        subsize[ic] = data->np_tot[i];
      } else if(i>dir) {
        startRecv[ic] = startSend[ic] = data->beg[i];
        subsize[ic] = data->np_int[i];
      } else {
        // That's the direction of exchange i==dir
        startRecv[ic] = (face == faceTop) ? 0 : data->end[i];
        startSend[ic] = (face == faceTop) ? data->beg[i] : data->end[i] - data->nghost[i];
        subsize[ic] = data->nghost[i];
      }
    }
    // We create two datatypes for these exchanges
    MPI_Datatype Send;
    MPI_Datatype Recv;
    MPI_SAFE_CALL(MPI_Type_create_subarray(3, size, subsize, startRecv,
                                            MPI_ORDER_C, realMPI, &Recv ));
    MPI_SAFE_CALL(MPI_Type_create_subarray(3, size, subsize, startSend,
                                            MPI_ORDER_C, realMPI, &Send ));
    MPI_SAFE_CALL(MPI_Type_commit(&Recv));
    MPI_SAFE_CALL(MPI_Type_commit(&Send));

    // now, stack together these datatypes for the variables we require
    std::vector<int> mapLengths(mapVars.size(), 1);    // Each block will be of size 1
    MPI_SAFE_CALL(MPI_Type_indexed(mapVars.size(), mapLengths.data(), mapVars.data(),
                                    Send, &typeVcSend[face] ));

    MPI_SAFE_CALL(MPI_Type_indexed(mapVars.size(), mapLengths.data(), mapVars.data(),
                                    Recv, &typeVcRecv[face] ));

    MPI_Type_commit(&typeVcSend[face]);
    MPI_Type_commit(&typeVcRecv[face]);

    // Free the remaining types
    MPI_Type_free(&Send);
    MPI_Type_free(&Recv);
  }
  #if MHD==YES
    // Array offsets (one additional element for each active dimension)
    int offset[3];
    offset[0] = IOFFSET;
    offset[1] = JOFFSET;
    offset[2] = KOFFSET;

    typeVsSend = std::vector<MPI_Datatype>(2);
    typeVsRecv = std::vector<MPI_Datatype>(2);
    // dir is the direction of exchange
    for(int face = 0 ; face <= 1 ; face++) {
      // We need to define one sub-array for each field component
      std::vector<MPI_Datatype> SendArray;
      std::vector<MPI_Datatype> RecvArray;
      std::vector<MPI_Aint> dispArray;
      MPI_Aint dispTotal = 0;
      for(int component = 0 ; component < DIMENSIONS ; component++) {
        MPI_Aint disp = 1;
        for(int i = 0 ; i < 3 ; i++) {
          int ic = 2-i; // because X1 is actually the last index in our arrays
          size[ic] = data->np_tot[i]+offset[i];
          disp = disp * size[ic];
          if(i<dir) {
            startRecv[ic] = startSend[ic] = 0;
            subsize[ic] = data->np_tot[i];
            if(i == component) subsize[ic]++;
          } else if(i>dir) {
            startRecv[ic] = startSend[ic] = data->beg[i];
            subsize[ic] = data->np_int[i];
            if(i == component) subsize[ic]++;
          } else {
            // That's the direction of exchange i==dir
            startRecv[ic] = (face == faceTop) ? 0 : data->end[i];
            startSend[ic] = (face == faceTop) ? data->beg[i] : data->end[i] - data->nghost[i];
            subsize[ic] = data->nghost[i];
          }
        }
        // We don't exchange the component that is normal to the direction dir
        if(component == dir) {
          dispTotal += disp;
          continue;
        }
        dispArray.push_back(dispTotal*sizeof(real));
        dispTotal += disp;
        SendArray.emplace_back();
        RecvArray.emplace_back();
        MPI_SAFE_CALL(MPI_Type_create_subarray(3, size, subsize, startRecv,
                                          MPI_ORDER_C, realMPI, &RecvArray.back() ));
        MPI_SAFE_CALL(MPI_Type_create_subarray(3, size, subsize, startSend,
                                          MPI_ORDER_C, realMPI, &SendArray.back() ));

        MPI_SAFE_CALL(MPI_Type_commit(&RecvArray.back()));
        MPI_SAFE_CALL(MPI_Type_commit(&SendArray.back()));
      }
      // SendArray and RecvArray now contains all of the arrays
      std::vector<int> arrayLength(RecvArray.size(), 1);    // Each block will be of size 1
      MPI_SAFE_CALL(
        MPI_Type_create_struct(RecvArray.size(), arrayLength.data(),
                                dispArray.data(), RecvArray.data(), &typeVsRecv[face]));
      MPI_SAFE_CALL(
        MPI_Type_create_struct(SendArray.size(), arrayLength.data(),
                                dispArray.data(), SendArray.data(), &typeVsSend[face]));

      MPI_Type_commit(&typeVsSend[face]);
      MPI_Type_commit(&typeVsRecv[face]);
      // Free the send/receive array per component on that face
      for(int i = 0 ; i < RecvArray.size() ; i++) {
        MPI_Type_free(&RecvArray[i]);
        MPI_Type_free(&SendArray[i]);
      }
    }
#endif // MHD
#endif // MPI
  idfx::popRegion();
}
