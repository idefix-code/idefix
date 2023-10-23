// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_CONSTRAINEDTRANSPORT_EMFEXCHANGE_HPP_
#define FLUID_CONSTRAINEDTRANSPORT_EMFEXCHANGE_HPP_

#include "constrainedTransport.hpp"
#include "fluid.hpp"
#include "dataBlock.hpp"

#ifdef WITH_MPI
template<typename Phys>
void ConstrainedTransport<Phys>::ExchangeAll() {
  if(data->mygrid->nproc[IDIR]>1) this->ExchangeX1();
  if(data->mygrid->nproc[JDIR]>1) this->ExchangeX2();
  if(data->mygrid->nproc[KDIR]>1) this->ExchangeX3();
}


// Exchange EMFs in X1
template<typename Phys>
void ConstrainedTransport<Phys>::ExchangeX1() {
  idfx::pushRegion("Emf::ExchangeX1");


  // Load  the buffers with data
  int ileft,iright,jbeg,jend,kbeg,kend;
  int ny;
  [[maybe_unused]] int nz;

  IdefixArray1D<real> BufferLeft=BufferSendX1[faceLeft];
  IdefixArray1D<real> BufferRight=BufferSendX1[faceRight];
  IdefixArray3D<real> ey=this->ey;
  IdefixArray3D<real> ez=this->ez;


  // If MPI Persistent, start receiving even before the buffers are filled
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];

  double tStart = MPI_Wtime();
  MPI_SAFE_CALL(MPI_Startall(2, recvRequestX1));
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;

  BoundaryType lbound = data->lbound[IDIR];
  BoundaryType rbound = data->rbound[IDIR];

  // Coordinates of the ghost region which needs to be transfered
  ileft   = data->beg[IDIR];
  iright = data->end[IDIR];
  jbeg   = data->beg[JDIR];
  jend   = data->end[JDIR];
  ny     = jend - jbeg;
  kbeg   = data->beg[KDIR];
  kend   = data->end[KDIR];
  nz     = kend - kbeg;

  idefix_for("LoadBufferX1Emfz",kbeg,kend,jbeg,jend+1,
    KOKKOS_LAMBDA (int k, int j) {
      BufferLeft( (j-jbeg) + (k-kbeg)*(ny+1) ) = ez(k,j,ileft);
      BufferRight( (j-jbeg) + (k-kbeg)*(ny+1) ) = ez(k,j,iright);
    }
  );
  #if DIMENSIONS == 3
  int Vsindex = (ny+1)*nz;

  idefix_for("LoadBufferX1Emfy",kbeg,kend+1,jbeg,jend,
    KOKKOS_LAMBDA (int k, int j) {
      BufferLeft( (j-jbeg) + (k-kbeg)*ny + Vsindex ) = ey(k,j,ileft);
      BufferRight( (j-jbeg) + (k-kbeg)*ny + Vsindex ) = ey(k,j,iright);
    }
  );
  #endif

  // Wait for completion before sending out everything
  Kokkos::fence();

  tStart = MPI_Wtime();
  MPI_SAFE_CALL(MPI_Startall(2, sendRequestX1));
  // Wait for buffers to be received

  MPI_Waitall(2,recvRequestX1,recvStatus);
  MPI_Waitall(2, sendRequestX1, sendStatus);
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;

  // Unpack
  BufferLeft=BufferRecvX1[faceLeft];
  BufferRight=BufferRecvX1[faceRight];

  // We average the edge emfs zones
  idefix_for("StoreBufferX1Emfz",kbeg,kend,jbeg,jend+1,
    KOKKOS_LAMBDA (int k, int j) {
      if(lbound == internal || lbound == periodic) {
        ez(k,j,ileft) = HALF_F*(
                        BufferLeft( (j-jbeg) + (k-kbeg)*(ny+1) ) + ez(k,j,ileft) );
      }
      if(rbound == internal || rbound == periodic) {
        ez(k,j,iright) = HALF_F*(
                        BufferRight( (j-jbeg) + (k-kbeg)*(ny+1) ) + ez(k,j,iright) );
      }
    });
  #if DIMENSIONS == 3
  Vsindex = (ny+1)*nz;
  idefix_for("StoreBufferX1Emfy",kbeg,kend+1,jbeg,jend,
    KOKKOS_LAMBDA (int k, int j) {
      if(lbound == internal || lbound == periodic) {
        ey(k,j,ileft) = HALF_F*(
                        BufferLeft( (j-jbeg) + (k-kbeg)*ny +Vsindex) + ey(k,j,ileft) );
      }
      if(rbound == internal || rbound == periodic) {
        ey(k,j,iright) = HALF_F*(
                        BufferRight( (j-jbeg) + (k-kbeg)*ny +Vsindex) + ey(k,j,iright) );
      }
    });
  #endif


  idfx::popRegion();
}

// Exchange EMFs in X2
template<typename Phys>
void ConstrainedTransport<Phys>::ExchangeX2() {
  idfx::pushRegion("Emf::ExchangeX2");

  // Load  the buffers with data
  int jleft,jright,ibeg,iend,kbeg,kend;
  int nx;
  [[maybe_unused]] int nz;
  IdefixArray1D<real> BufferLeft=BufferSendX2[faceLeft];
  IdefixArray1D<real> BufferRight=BufferSendX2[faceRight];
  IdefixArray3D<real> ex=this->ex;
  IdefixArray3D<real> ez=this->ez;

  // If MPI Persistent, start receiving even before the buffers are filled
  double tStart = MPI_Wtime();
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];
  MPI_SAFE_CALL(MPI_Startall(2, recvRequestX2));
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;

  BoundaryType lbound = data->lbound[JDIR];
  BoundaryType rbound = data->rbound[JDIR];

  // Coordinates of the ghost region which needs to be transfered
  ibeg   = data->beg[IDIR];
  iend   = data->end[IDIR];
  nx     = iend - ibeg;
  jleft   = data->beg[JDIR];
  jright = data->end[JDIR];
  kbeg   = data->beg[KDIR];
  kend   = data->end[KDIR];
  nz     = kend - kbeg;

  idefix_for("LoadBufferX2Emfz",kbeg,kend,ibeg,iend+1,
    KOKKOS_LAMBDA (int k, int i) {
      BufferLeft( (i-ibeg) + (k-kbeg)*(nx+1) ) = ez(k,jleft,i);
      BufferRight( (i-ibeg) + (k-kbeg)*(nx+1) ) = ez(k,jright,i);
    }
  );
  #if DIMENSIONS == 3
  int Vsindex = (nx+1)*nz;

  idefix_for("LoadBufferX1Emfx",kbeg,kend+1,ibeg,iend,
    KOKKOS_LAMBDA (int k, int i) {
      BufferLeft( (i-ibeg) + (k-kbeg)*nx + Vsindex ) = ex(k,jleft,i);
      BufferRight( (i-ibeg) + (k-kbeg)*nx + Vsindex ) = ex(k,jright,i);
    }
  );
  #endif

  // Wait for completion before sending out everything
  Kokkos::fence();

  tStart = MPI_Wtime();
  MPI_SAFE_CALL(MPI_Startall(2, sendRequestX2));
  // Wait for buffers to be received
  MPI_Waitall(2,recvRequestX2,recvStatus);
  MPI_Waitall(2, sendRequestX2, sendStatus);
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;

  // Unpack
  BufferLeft=BufferRecvX2[faceLeft];
  BufferRight=BufferRecvX2[faceRight];

  // We average the edge emfs zones
  idefix_for("StoreBufferX2Emfz",kbeg,kend,ibeg,iend+1,
    KOKKOS_LAMBDA (int k, int i) {
      if(lbound == internal || lbound == periodic) {
        ez(k,jleft,i) = HALF_F*(
                        BufferLeft( (i-ibeg) + (k-kbeg)*(nx+1) ) + ez(k,jleft,i) );
      }
      if(rbound == internal || rbound == periodic) {
        ez(k,jright,i) = HALF_F*(
                        BufferRight( (i-ibeg) + (k-kbeg)*(nx+1) ) + ez(k,jright,i) );
      }
    });
  #if DIMENSIONS == 3
  Vsindex = (nx+1)*nz;
  idefix_for("StoreBufferX1Emfy",kbeg,kend+1,ibeg,iend,
    KOKKOS_LAMBDA (int k, int i) {
      if(lbound == internal || lbound == periodic) {
        ex(k,jleft,i) = HALF_F*(
                        BufferLeft( (i-ibeg) + (k-kbeg)*nx +Vsindex) + ex(k,jleft,i) );
      }
      if(rbound == internal || rbound == periodic) {
        ex(k,jright,i) = HALF_F*(
                        BufferRight( (i-ibeg) + (k-kbeg)*nx +Vsindex) + ex(k,jright,i) );
      }
    });
  #endif


  idfx::popRegion();
}

// Exchange EMFs in X3
template<typename Phys>
void ConstrainedTransport<Phys>::ExchangeX3() {
  idfx::pushRegion("Emf::ExchangeX3");


  // Load  the buffers with data
  int kleft,kright,ibeg,iend,jbeg,jend;
  int nx,ny;
  IdefixArray1D<real> BufferLeft=BufferSendX3[faceLeft];
  IdefixArray1D<real> BufferRight=BufferSendX3[faceRight];
  IdefixArray3D<real> ex=this->ex;
  IdefixArray3D<real> ey=this->ey;

  int Vsindex = 0;


  // If MPI Persistent, start receiving even before the buffers are filled
  double tStart = MPI_Wtime();
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];
  MPI_SAFE_CALL(MPI_Startall(2, recvRequestX3));
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;

  BoundaryType lbound = data->lbound[KDIR];
  BoundaryType rbound = data->rbound[KDIR];

  // Coordinates of the ghost region which needs to be transfered
  ibeg   = data->beg[IDIR];
  iend   = data->end[IDIR];
  nx     = iend - ibeg;

  jbeg   = data->beg[JDIR];
  jend   = data->end[JDIR];
  ny     = jend - jbeg;

  kleft   = data->beg[KDIR];
  kright = data->end[KDIR];

  idefix_for("LoadBufferX3Emfx",jbeg,jend+1,ibeg,iend,
    KOKKOS_LAMBDA (int j, int i) {
      BufferLeft( (i-ibeg) + (j-jbeg)*nx ) = ex(kleft,j,i);
      BufferRight( (i-ibeg) + (j-jbeg)*nx ) = ex(kright,j,i);
    }
  );
  Vsindex = nx*(ny+1);

  idefix_for("LoadBufferX3Emfy",jbeg,jend,ibeg,iend+1,
    KOKKOS_LAMBDA (int j, int i) {
      BufferLeft( (i-ibeg) + (j-jbeg)*(nx+1) + Vsindex ) = ey(kleft,j,i);
      BufferRight( (i-ibeg) + (j-jbeg)*(nx+1) + Vsindex ) = ey(kright,j,i);
    }
  );

  // Wait for completion before sending out everything
  Kokkos::fence();

  tStart = MPI_Wtime();
  MPI_SAFE_CALL(MPI_Startall(2, sendRequestX3));
  // Wait for buffers to be received
  MPI_Waitall(2,recvRequestX3,recvStatus);
  MPI_Waitall(2, sendRequestX3, sendStatus);
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;

  // Unpack
  BufferLeft=BufferRecvX3[faceLeft];
  BufferRight=BufferRecvX3[faceRight];

  // We average the edge emfs zones
  idefix_for("StoreBufferX3Emfx",jbeg,jend+1,ibeg,iend,
    KOKKOS_LAMBDA (int j, int i) {
      if(lbound == internal || lbound == periodic) {
        ex(kleft,j,i) = HALF_F*(
                        BufferLeft( (i-ibeg) + (j-jbeg)*nx ) + ex(kleft,j,i) );
      }
      if(rbound == internal || rbound == periodic) {
        ex(kright,j,i) = HALF_F*(
                        BufferRight( (i-ibeg) + (j-jbeg)*nx ) + ex(kright,j,i) );
      }
    });

  Vsindex = nx*(ny+1);
  idefix_for("StoreBufferX3Emfy",jbeg,jend,ibeg,iend+1,
    KOKKOS_LAMBDA (int j, int i) {
      if(lbound == internal || lbound == periodic) {
        ey(kleft,j,i) = HALF_F*(
                        BufferLeft( (i-ibeg) + (j-jbeg)*(nx+1) + Vsindex ) + ey(kleft,j,i) );
      }
      if(rbound == internal || rbound == periodic) {
        ey(kright,j,i) = HALF_F*(
                        BufferRight( (i-ibeg) + (j-jbeg)*(nx+1) + Vsindex ) + ey(kright,j,i) );
      }
    });


  idfx::popRegion();
}

#endif
#endif // FLUID_CONSTRAINEDTRANSPORT_EMFEXCHANGE_HPP_
