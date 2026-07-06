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
  int ileft,iright,jbeg,jend;

  Buffer BufferLeft=BufferSendX1[faceLeft];
  Buffer BufferRight=BufferSendX1[faceRight];


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
  BufferRight.Pack(ez, sendBoxEz);
  #if DIMENSIONS == 3

  //extend by one the end on kdir && take the ghost on i
  BoundingBox sendBoxEy = baseBox;
  sendBoxEy[KDIR][1] += 1;
  sendBoxEy[IDIR][0] = iright;
  sendBoxEy[IDIR][1] = iright + 1;
  BufferRight.Pack(ez, sendBoxEy);
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

  // Erase the emf with the one coming from the left process
  //extend by one the end on jdir && take the ghost revc zone on i
  BoundingBox recvBoxEz = baseBox;
  recvBoxEz[JDIR][1] += 1;
  recvBoxEz[IDIR][0] = ileft;
  recvBoxEz[IDIR][1] = ileft + 1;
  if(lbound == internal || lbound == periodic)
    BufferLeft.Unpack(ez, recvBoxEz);

  #if DIMENSIONS == 3
  //extend by one the end on kdir && take the ghost revc zone on i
  BoundingBox recvBoxEy = baseBox;
  recvBoxEy[KDIR][1] += 1;
  recvBoxEy[IDIR][0] = ileft;
  recvBoxEy[IDIR][1] = ileft + 1;
  if(lbound == internal || lbound == periodic)
    BufferLeft.Unpack(ey, recvBoxEy);
  #endif


  idfx::popRegion();
}

// Exchange EMFs in X2
template<typename Phys>
void ConstrainedTransport<Phys>::ExchangeX2(IdefixArray3D<real> ex, IdefixArray3D<real> ez) {
  idfx::pushRegion("Emf::ExchangeX2");

  // Load  the buffers with data
  int jleft,jright;
  Buffer BufferLeft=BufferSendX2[faceLeft];
  Buffer BufferRight=BufferSendX2[faceRight];

  // If MPI Persistent, start receiving even before the buffers are filled
  double tStart = MPI_Wtime();
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];
  MPI_SAFE_CALL(MPI_Startall(2, recvRequestX2));
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;

  BoundaryType lbound = data->lbound[JDIR];
  BoundaryType rbound = data->rbound[JDIR];

  // Coordinates of the ghost region which needs to be transfered
  jleft   = data->beg[JDIR];
  jright = data->end[JDIR];

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
  BufferRight.Pack(ez, sendBoxEz);
  #if DIMENSIONS == 3

  //extend by one the end on idir && take the ghost on j
  BoundingBox sendBoxEx = baseBox;
  sendBoxEx[KDIR][1] += 1;
  sendBoxEx[JDIR][0] = jright;
  sendBoxEx[JDIR][1] = jright + 1;
  BufferRight.Pack(ex, sendBoxEx);
  #endif

  // Wait for completion before sending out everything
  Kokkos::fence();

  tStart = MPI_Wtime();
  MPI_SAFE_CALL(MPI_Startall(2, sendRequestX2));
  // Wait for buffers to be received
  MPI_Waitall(2, recvRequestX2, recvStatus);
  MPI_Waitall(2, sendRequestX2, sendStatus);
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;

  // Unpack
  BufferLeft=BufferRecvX2[faceLeft];
  BufferRight=BufferRecvX2[faceRight];

  // We average the edge emfs zones
  //extend by one the end on idir && take the ghost revc zone on j
  BoundingBox recvBoxEz = baseBox;
  recvBoxEz[IDIR][1] += 1;
  recvBoxEz[JDIR][0] = jleft;
  recvBoxEz[JDIR][1] = jleft + 1;
  if(lbound == internal || lbound == periodic)
    BufferLeft.Unpack(ez, recvBoxEz);
  #if DIMENSIONS == 3
  //extend by one the end on kdir && take the ghost revc zone on j
  BoundingBox recvBoxEx = baseBox;
  recvBoxEx[KDIR][1] += 1;
  recvBoxEx[JDIR][0] = jleft;
  recvBoxEx[JDIR][1] = jleft + 1;
  if(lbound == internal || lbound == periodic)
    BufferLeft.Unpack(ex, recvBoxEx);
  #endif


  idfx::popRegion();
}

// Exchange EMFs in X3
template<typename Phys>
void ConstrainedTransport<Phys>::ExchangeX3(IdefixArray3D<real> ex, IdefixArray3D<real> ey) {
  idfx::pushRegion("Emf::ExchangeX3");


  // Load  the buffers with data
  int kleft,kright,jbeg,jend;
  Buffer BufferLeft=BufferSendX3[faceLeft];
  Buffer BufferRight=BufferSendX3[faceRight];

  // If MPI Persistent, start receiving even before the buffers are filled
  double tStart = MPI_Wtime();
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];
  MPI_SAFE_CALL(MPI_Startall(2, recvRequestX3));
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;

  BoundaryType lbound = data->lbound[KDIR];
  BoundaryType rbound = data->rbound[KDIR];

  // Coordinates of the ghost region which needs to be transfered
  jbeg   = data->beg[JDIR];
  jend   = data->end[JDIR];

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
  BufferRight.Pack(ez, sendBoxEz);

  //extend by one the end on idir && take the ghost on k
  BoundingBox sendBoxEy = baseBox;
  sendBoxEy[IDIR][1] += 1;
  sendBoxEy[KDIR][0] = kright;
  sendBoxEy[KDIR][1] = kright + 1;
  BufferRight.Pack(ey, sendBoxEy);

  // Wait for completion before sending out everything
  Kokkos::fence();

  tStart = MPI_Wtime();
  MPI_SAFE_CALL(MPI_Startall(2, sendRequestX3));
  // Wait for buffers to be received
  MPI_Waitall(2, recvRequestX3, recvStatus);
  MPI_Waitall(2, sendRequestX3, sendStatus);
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;

  // Unpack
  BufferLeft=BufferRecvX3[faceLeft];
  BufferRight=BufferRecvX3[faceRight];

  // We average the edge emfs zones
  //extend by one the end on jdir && take the ghost revc zone on k
  BoundingBox recvBoxEx = baseBox;
  recvBoxEx[JDIR][1] += 1;
  recvBoxEx[KDIR][0] = kleft;
  recvBoxEx[KDIR][1] = kleft + 1;
  if(lbound == internal || lbound == periodic)
    BufferLeft.Unpack(ex, recvBoxEx);

  //extend by one the end on idir && take the ghost revc zone on k
  BoundingBox recvBoxEy = baseBox;
  recvBoxEy[IDIR][1] += 1;
  recvBoxEy[KDIR][0] = kleft;
  recvBoxEy[KDIR][1] = kleft + 1;
  if(lbound == internal || lbound == periodic)
    BufferLeft.Unpack(ey, recvBoxEy);
  idfx::popRegion();
}

#endif
#endif // FLUID_CONSTRAINEDTRANSPORT_EMFEXCHANGE_HPP_
