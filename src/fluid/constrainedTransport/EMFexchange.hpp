// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_CONSTRAINEDTRANSPORT_EMFEXCHANGE_HPP_
#define FLUID_CONSTRAINEDTRANSPORT_EMFEXCHANGE_HPP_

#include <cassert>

#include "constrainedTransport.hpp"
#include "fluid.hpp"
#include "dataBlock.hpp"

#ifdef WITH_MPI
template<typename Phys>
void ConstrainedTransport<Phys>::ExchangeAll(IdefixArray3D<real> ex,
                                             IdefixArray3D<real> ey,
                                             IdefixArray3D<real> ez) {
  if(data->mygrid->nproc[IDIR]>1) this->ExchangeX1(ey,ez);
  if(data->mygrid->nproc[JDIR]>1) this->ExchangeX2(ex,ez);
  if(data->mygrid->nproc[KDIR]>1) this->ExchangeX3(ex,ey);
}


// Exchange EMFs in X1
template<typename Phys>
void ConstrainedTransport<Phys>::ExchangeX1(IdefixArray3D<real> ey, IdefixArray3D<real> ez) {
  idfx::pushRegion("Emf::ExchangeX1");


  // Load  the buffers with data
  int ileft,iright,jbeg,jend,kbeg,kend;
  int ny;
  [[maybe_unused]] int nz;

  IdefixArray1D<real> BufferLeft=BufferSendX1[faceLeft];
  IdefixArray1D<real> BufferRight=BufferSendX1[faceRight];
  Buffer BufferLeftNew=BufferSendX1New[faceLeft];
  Buffer BufferRightNew=BufferSendX1New[faceRight];


  // If MPI Persistent, start receiving even before the buffers are filled
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];
  MPI_Status sendStatusNew[2];
  MPI_Status recvStatusNew[2];

  double tStart = MPI_Wtime();
  MPI_SAFE_CALL(MPI_Startall(2, recvRequestX1));
  MPI_SAFE_CALL(MPI_Startall(2, recvRequestX1New));
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

  // Create the base box to be patched by the ops
  BoundingBox baseBox;
  baseBox[IDIR][0] = data->beg[IDIR];
  baseBox[IDIR][1] = data->end[IDIR];
  baseBox[JDIR][0] = data->beg[JDIR];
  baseBox[JDIR][1] = data->end[JDIR];
  baseBox[KDIR][0] = data->beg[KDIR];
  baseBox[KDIR][1] = data->end[KDIR];

  //extend by one the end on jdir && take the ghost on i
  BoundingBox sendBoxEz = baseBox;
  sendBoxEz[JDIR][1] += 1;
  sendBoxEz[IDIR][0] = iright;
  sendBoxEz[IDIR][1] = iright + 1;
  BufferRightNew.Pack(ez, sendBoxEz);
  idefix_for("LoadBufferX1Emfz",kbeg,kend,jbeg,jend+1,
    KOKKOS_LAMBDA (int k, int j) {
      BufferRight( (j-jbeg) + (k-kbeg)*(ny+1) ) = ez(k,j,iright);
      assert(BufferRightNew.getArray()((j-jbeg) + (k-kbeg)*(ny+1)) == ez(k,j,iright));
    }
  );
  #if DIMENSIONS == 3
  int Vsindex = (ny+1)*nz;

  //extend by one the end on kdir && take the ghost on i
  BoundingBox sendBoxEy = baseBox;
  sendBoxEy[KDIR][1] += 1;
  sendBoxEy[IDIR][0] = iright;
  sendBoxEy[IDIR][1] = iright + 1;
  BufferRightNew.Pack(ez, sendBoxEy);
  idefix_for("LoadBufferX1Emfy",kbeg,kend+1,jbeg,jend,
    KOKKOS_LAMBDA (int k, int j) {
      BufferRight( (j-jbeg) + (k-kbeg)*ny + Vsindex ) = ey(k,j,iright);
      assert(BufferRightNew.getArray()( (j-jbeg) + (k-kbeg)*ny + Vsindex ) == ey(k,j,iright));
    }
  );
  #endif

  assert(BufferRight == BufferRightNew.getArray());

  // Wait for completion before sending out everything
  Kokkos::fence();

  tStart = MPI_Wtime();
  MPI_SAFE_CALL(MPI_Startall(2, sendRequestX1));
  MPI_SAFE_CALL(MPI_Startall(2, sendRequestX1New));
  // Wait for buffers to be received

  MPI_Waitall(2,recvRequestX1,recvStatus);
  MPI_Waitall(2,recvRequestX1New,recvStatusNew);
  MPI_Waitall(2, sendRequestX1, sendStatus);
  MPI_Waitall(2, sendRequestX1New, sendStatusNew);
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;

  // Unpack
  BufferLeft=BufferRecvX1[faceLeft];
  BufferRight=BufferRecvX1[faceRight];
  BufferLeftNew=BufferRecvX1New[faceLeft];
  BufferRightNew=BufferRecvX1New[faceRight];



  assert(BufferLeft == BufferLeftNew.getArray());

  // Erase the emf with the one coming from the left process
  //extend by one the end on jdir && take the ghost revc zone on i
  BoundingBox recvBoxEz = baseBox;
  recvBoxEz[JDIR][1] += 1;
  recvBoxEz[IDIR][0] = ileft;
  recvBoxEz[IDIR][1] = ileft + 1;
  if(lbound == internal || lbound == periodic)
    BufferLeftNew.Unpack(ez, recvBoxEz);
  idefix_for("StoreBufferX1Emfz",kbeg,kend,jbeg,jend+1,
    KOKKOS_LAMBDA (int k, int j) {
      if(lbound == internal || lbound == periodic) {//////////////<<<<<<<<<<<<<<<<<<TODO
        assert(ez(k,j,ileft) == BufferLeft( (j-jbeg) + (k-kbeg)*(ny+1)));
        ez(k,j,ileft) = BufferLeft( (j-jbeg) + (k-kbeg)*(ny+1));
      }
    });

  #if DIMENSIONS == 3
  Vsindex = (ny+1)*nz;
  //extend by one the end on kdir && take the ghost revc zone on i
  BoundingBox recvBoxEy = baseBox;
  recvBoxEy[KDIR][1] += 1;
  recvBoxEy[IDIR][0] = ileft;
  recvBoxEy[IDIR][1] = ileft + 1;
  if(lbound == internal || lbound == periodic)
    BufferLeftNew.Unpack(ey, recvBoxEy);
  idefix_for("StoreBufferX1Emfy",kbeg,kend+1,jbeg,jend,
    KOKKOS_LAMBDA (int k, int j) {
      if(lbound == internal || lbound == periodic) {//////////////<<<<<<<<<<<<<<<<<<TODO
        assert(ey(k,j,ileft) == BufferLeft( (j-jbeg) + (k-kbeg)*ny +Vsindex));
        ey(k,j,ileft) = BufferLeft( (j-jbeg) + (k-kbeg)*ny +Vsindex);
      }
    });
  #endif


  idfx::popRegion();
}

// Exchange EMFs in X2
template<typename Phys>
void ConstrainedTransport<Phys>::ExchangeX2(IdefixArray3D<real> ex, IdefixArray3D<real> ez) {
  idfx::pushRegion("Emf::ExchangeX2");

  // Load  the buffers with data
  int jleft,jright,ibeg,iend,kbeg,kend;
  int nx;
  [[maybe_unused]] int nz;
  IdefixArray1D<real> BufferLeft=BufferSendX2[faceLeft];
  IdefixArray1D<real> BufferRight=BufferSendX2[faceRight];
  Buffer BufferLeftNew=BufferSendX2New[faceLeft];
  Buffer BufferRightNew=BufferSendX2New[faceRight];

  // If MPI Persistent, start receiving even before the buffers are filled
  double tStart = MPI_Wtime();
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];
  MPI_Status sendStatusNew[2];
  MPI_Status recvStatusNew[2];
  MPI_SAFE_CALL(MPI_Startall(2, recvRequestX2));
  MPI_SAFE_CALL(MPI_Startall(2, recvRequestX2New));
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

  // Create the base box to be patched by the ops
  BoundingBox baseBox;
  baseBox[IDIR][0] = data->beg[IDIR];
  baseBox[IDIR][1] = data->end[IDIR];
  baseBox[JDIR][0] = data->beg[JDIR];
  baseBox[JDIR][1] = data->end[JDIR];
  baseBox[KDIR][0] = data->beg[KDIR];
  baseBox[KDIR][1] = data->end[KDIR];

  //extend by one the end on idir && take the ghost on j
  BoundingBox sendBoxEz = baseBox;
  sendBoxEz[IDIR][1] += 1;
  sendBoxEz[JDIR][0] = jright;
  sendBoxEz[JDIR][1] = jright + 1;
  BufferRightNew.Pack(ez, sendBoxEz);
  idefix_for("LoadBufferX2Emfz",kbeg,kend,ibeg,iend+1,
    KOKKOS_LAMBDA (int k, int i) {
      BufferRight( (i-ibeg) + (k-kbeg)*(nx+1) ) = ez(k,jright,i);
      assert(BufferRightNew.getArray()( (i-ibeg) + (k-kbeg)*(nx+1) ) == ez(k,jright,i));
    }
  );
  #if DIMENSIONS == 3
  int Vsindex = (nx+1)*nz;

  //extend by one the end on idir && take the ghost on j
  BoundingBox sendBoxEx = baseBox;
  sendBoxEx[KDIR][1] += 1;
  sendBoxEx[JDIR][0] = jright;
  sendBoxEx[JDIR][1] = jright + 1;
  BufferRightNew.Pack(ex, sendBoxEx);
  idefix_for("LoadBufferX1Emfx",kbeg,kend+1,ibeg,iend,
    KOKKOS_LAMBDA (int k, int i) {
      BufferRight( (i-ibeg) + (k-kbeg)*nx + Vsindex ) = ex(k,jright,i);
      assert(BufferRightNew.getArray()( (i-ibeg) + (k-kbeg)*nx + Vsindex ) == ex(k,jright,i));
    }
  );
  #endif

  assert(BufferRightNew.getArray() == BufferRight);

  // Wait for completion before sending out everything
  Kokkos::fence();

  tStart = MPI_Wtime();
  MPI_SAFE_CALL(MPI_Startall(2, sendRequestX2));
  MPI_SAFE_CALL(MPI_Startall(2, sendRequestX2New));
  // Wait for buffers to be received
  MPI_Waitall(2,recvRequestX2,recvStatus);
  MPI_Waitall(2,recvRequestX2New,recvStatus);
  MPI_Waitall(2, sendRequestX2, sendStatus);
  MPI_Waitall(2, sendRequestX2New, sendStatus);
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;

  // Unpack
  BufferLeft=BufferRecvX2[faceLeft];
  BufferRight=BufferRecvX2[faceRight];
  BufferLeftNew=BufferRecvX2New[faceLeft];
  BufferRightNew=BufferRecvX2New[faceRight];

  assert(BufferLeftNew.getArray() == BufferLeft);

  // We average the edge emfs zones
  //extend by one the end on idir && take the ghost revc zone on j
  BoundingBox recvBoxEz = baseBox;
  recvBoxEz[IDIR][1] += 1;
  recvBoxEz[JDIR][0] = jleft;
  recvBoxEz[JDIR][1] = jleft + 1;
  if(lbound == internal || lbound == periodic)
    BufferLeftNew.Unpack(ez, recvBoxEz);
  idefix_for("StoreBufferX2Emfz",kbeg,kend,ibeg,iend+1,
    KOKKOS_LAMBDA (int k, int i) {
      if(lbound == internal || lbound == periodic) {//////////////<<<<<<<<<<<<<<<<<<TODO
        assert(ez(k,jleft,i) == BufferLeft( (i-ibeg) + (k-kbeg)*(nx+1) ));
        ez(k,jleft,i) = BufferLeft( (i-ibeg) + (k-kbeg)*(nx+1) );
      }
    });
  #if DIMENSIONS == 3
  Vsindex = (nx+1)*nz;
  //extend by one the end on kdir && take the ghost revc zone on j
  BoundingBox recvBoxEx = baseBox;
  recvBoxEx[KDIR][1] += 1;
  recvBoxEx[JDIR][0] = jleft;
  recvBoxEx[JDIR][1] = jleft + 1;
  if(lbound == internal || lbound == periodic)
    BufferLeftNew.Unpack(ex, recvBoxEx);
  idefix_for("StoreBufferX1Emfy",kbeg,kend+1,ibeg,iend,
    KOKKOS_LAMBDA (int k, int i) {
      if(lbound == internal || lbound == periodic) {//////////////<<<<<<<<<<<<<<<<<<TODO
        assert(ex(k,jleft,i) == BufferLeft( (i-ibeg) + (k-kbeg)*nx +Vsindex));
        ex(k,jleft,i) = BufferLeft( (i-ibeg) + (k-kbeg)*nx +Vsindex);
      }
    });
  #endif


  idfx::popRegion();
}

// Exchange EMFs in X3
template<typename Phys>
void ConstrainedTransport<Phys>::ExchangeX3(IdefixArray3D<real> ex, IdefixArray3D<real> ey) {
  idfx::pushRegion("Emf::ExchangeX3");


  // Load  the buffers with data
  int kleft,kright,ibeg,iend,jbeg,jend;
  int nx,ny;
  IdefixArray1D<real> BufferLeft=BufferSendX3[faceLeft];
  IdefixArray1D<real> BufferRight=BufferSendX3[faceRight];
  Buffer BufferLeftNew=BufferSendX3New[faceLeft];
  Buffer BufferRightNew=BufferSendX3New[faceRight];

  int Vsindex = 0;


  // If MPI Persistent, start receiving even before the buffers are filled
  double tStart = MPI_Wtime();
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];
  MPI_Status sendStatusNew[2];
  MPI_Status recvStatusNew[2];
  MPI_SAFE_CALL(MPI_Startall(2, recvRequestX3));
  MPI_SAFE_CALL(MPI_Startall(2, recvRequestX3New));
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

  // Create the base box to be patched by the ops
  BoundingBox baseBox;
  baseBox[IDIR][0] = data->beg[IDIR];
  baseBox[IDIR][1] = data->end[IDIR];
  baseBox[JDIR][0] = data->beg[JDIR];
  baseBox[JDIR][1] = data->end[JDIR];
  baseBox[KDIR][0] = data->beg[KDIR];
  baseBox[KDIR][1] = data->end[KDIR];

  //extend by one the end on jdir && take the ghost on k
  BoundingBox sendBoxEz = baseBox;
  sendBoxEz[JDIR][1] += 1;
  sendBoxEz[KDIR][0] = kright;
  sendBoxEz[KDIR][1] = kright + 1;
  BufferRightNew.Pack(ez, sendBoxEz);
  idefix_for("LoadBufferX3Emfx",jbeg,jend+1,ibeg,iend,
    KOKKOS_LAMBDA (int j, int i) {
      BufferRight( (i-ibeg) + (j-jbeg)*nx ) = ex(kright,j,i);
      assert(BufferRightNew.getArray()( (i-ibeg) + (j-jbeg)*nx ) == ex(kright,j,i));
    }
  );
  Vsindex = nx*(ny+1);

  //extend by one the end on idir && take the ghost on k
  BoundingBox sendBoxEy = baseBox;
  sendBoxEy[IDIR][1] += 1;
  sendBoxEy[KDIR][0] = kright;
  sendBoxEy[KDIR][1] = kright + 1;
  BufferRightNew.Pack(ey, sendBoxEy);
  idefix_for("LoadBufferX3Emfy",jbeg,jend,ibeg,iend+1,
    KOKKOS_LAMBDA (int j, int i) {
      BufferRight( (i-ibeg) + (j-jbeg)*(nx+1) + Vsindex ) = ey(kright,j,i);
      assert(BufferRightNew.getArray()( (i-ibeg) + (j-jbeg)*(nx+1) + Vsindex ) == ey(kright,j,i));
    }
  );

  assert(BufferRightNew.getArray() == BufferRight);

  // Wait for completion before sending out everything
  Kokkos::fence();

  tStart = MPI_Wtime();
  MPI_SAFE_CALL(MPI_Startall(2, sendRequestX3));
  MPI_SAFE_CALL(MPI_Startall(2, sendRequestX3New));
  // Wait for buffers to be received
  MPI_Waitall(2,recvRequestX3,recvStatus);
  MPI_Waitall(2,recvRequestX3New,recvStatusNew);
  MPI_Waitall(2, sendRequestX3, sendStatus);
  MPI_Waitall(2, sendRequestX3New, sendStatusNew);
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;

  // Unpack
  BufferLeft=BufferRecvX3[faceLeft];
  BufferRight=BufferRecvX3[faceRight];
  BufferLeftNew=BufferRecvX3New[faceLeft];
  BufferRightNew=BufferRecvX3New[faceRight];

  assert(BufferLeftNew.getArray() == BufferLeft);

  // We average the edge emfs zones
  //extend by one the end on jdir && take the ghost revc zone on k
  BoundingBox recvBoxEx = baseBox;
  recvBoxEx[JDIR][1] += 1;
  recvBoxEx[KDIR][0] = kleft;
  recvBoxEx[KDIR][1] = kleft + 1;
  if(lbound == internal || lbound == periodic)
    BufferLeftNew.Unpack(ex, recvBoxEx);
  idefix_for("StoreBufferX3Emfx",jbeg,jend+1,ibeg,iend,
    KOKKOS_LAMBDA (int j, int i) {
      if(lbound == internal || lbound == periodic) {//////////////<<<<<<<<<<<<<<<<<<TODO
        assert(ex(kleft,j,i) == BufferLeft( (i-ibeg) + (j-jbeg)*nx ));
        ex(kleft,j,i) = BufferLeft( (i-ibeg) + (j-jbeg)*nx );
      }
    });

  Vsindex = nx*(ny+1);
  //extend by one the end on idir && take the ghost revc zone on k
  BoundingBox recvBoxEy = baseBox;
  recvBoxEy[IDIR][1] += 1;
  recvBoxEy[KDIR][0] = kleft;
  recvBoxEy[KDIR][1] = kleft + 1;
  if(lbound == internal || lbound == periodic)
    BufferLeftNew.Unpack(ey, recvBoxEy);
  idefix_for("StoreBufferX3Emfy",jbeg,jend,ibeg,iend+1,
    KOKKOS_LAMBDA (int j, int i) {
      if(lbound == internal || lbound == periodic) {//////////////<<<<<<<<<<<<<<<<<<TODO
        assert(ey(kleft,j,i) == BufferLeft( (i-ibeg) + (j-jbeg)*(nx+1) + Vsindex ));
        ey(kleft,j,i) = BufferLeft( (i-ibeg) + (j-jbeg)*(nx+1) + Vsindex );
      }
    });
  idfx::popRegion();
}

#endif
#endif // FLUID_CONSTRAINEDTRANSPORT_EMFEXCHANGE_HPP_
