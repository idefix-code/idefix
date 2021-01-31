// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "../../idefix.hpp"
#include "dataBlock.hpp"
#include "mpi.hpp"

//#define MPI_NON_BLOCKING
#define MPI_PERSISTENT

// MPI Routines exchange
void Mpi::ExchangeAll() {
  IDEFIX_ERROR("Not Implemented");
}

void Mpi::InitFromDataBlock(DataBlock *datain) {
  this->data = datain;
  this->mygrid = datain->mygrid;
  this->timer.reset();
  ////////////////////////////////////////////////////////////////////////////
  // Init variable mappers
  // The variable mapper allows associates on number to a variable number in Vc
  // This is required since we skip some of the variables in Vc
  #if MHD == YES
  mapNVars = NVAR - DIMENSIONS; // We don't send magnetic field components which are in Vs
  #else
  mapNVars = NVAR;
  #endif

  mapVars = IdefixArray1D<int>("mapVars",mapNVars);

  // Create a host mirror
  IdefixArray1D<int>::HostMirror mapVarsHost = Kokkos::create_mirror_view(mapVars);
  // Init the host mirror
  int ntarget = 0;
  for(int n = 0 ; n < mapNVars ; n++) {
    mapVarsHost(n) = ntarget;
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
  // Synchronize the mapVars
  Kokkos::deep_copy(mapVars,mapVarsHost);

  /////////////////////////////////////////////////////////////////////////////
  // Init exchange datasets
  bufferSizeX1 = 0;
  bufferSizeX2 = 0;
  bufferSizeX3 = 0;

  // Number of cells in X1 boundary condition:
  bufferSizeX1 = data->nghost[IDIR] * data->np_int[JDIR] * data->np_int[KDIR] * mapNVars;

#if MHD == YES
  // BX1s
  bufferSizeX1 += data->np_int[JDIR] * data->np_int[KDIR];

  #if DIMENSIONS>=2
  bufferSizeX1 += data->nghost[IDIR] * (data->np_int[JDIR]+1) * data->np_int[KDIR];
  #endif

  #if DIMENSIONS==3
  bufferSizeX1 += data->nghost[IDIR] * data->np_int[JDIR] * (data->np_int[KDIR]+1);
  #endif  // DIMENSIONS
#endif  // MHD

  BufferRecvX1[faceLeft ] = IdefixArray1D<real>("BufferRecvX1Left", bufferSizeX1);
  BufferRecvX1[faceRight] = IdefixArray1D<real>("BufferRecvX1Right",bufferSizeX1);
  BufferSendX1[faceLeft ] = IdefixArray1D<real>("BufferSendX1Left", bufferSizeX1);
  BufferSendX1[faceRight] = IdefixArray1D<real>("BufferSendX1Right",bufferSizeX1);

  // Number of cells in X2 boundary condition (only required when problem >2D):
#if DIMENSIONS >= 2
  bufferSizeX2 = data->np_tot[IDIR] * data->nghost[JDIR] * data->np_int[KDIR] * mapNVars;
  #if MHD == YES
  // BX1s
  bufferSizeX2 += (data->np_tot[IDIR]+1) * data->nghost[JDIR] * data->np_int[KDIR];
  // BX2s
  bufferSizeX2 += data->np_tot[IDIR] * data->np_int[KDIR];
    #if DIMENSIONS==3
  bufferSizeX2 += data->np_tot[IDIR] * data->nghost[JDIR] * (data->np_int[KDIR]+1);
    #endif  // DIMENSIONS
  #endif  // MHD

  BufferRecvX2[faceLeft ] = IdefixArray1D<real>("BufferRecvX2Left", bufferSizeX2);
  BufferRecvX2[faceRight] = IdefixArray1D<real>("BufferRecvX2Right",bufferSizeX2);
  BufferSendX2[faceLeft ] = IdefixArray1D<real>("BufferSendX2Left", bufferSizeX2);
  BufferSendX2[faceRight] = IdefixArray1D<real>("BufferSendX2Right",bufferSizeX2);

#endif
// Number of cells in X3 boundary condition (only required when problem is 3D):
#if DIMENSIONS ==3
  bufferSizeX3 = data->np_tot[IDIR] * data->np_tot[JDIR] * data->nghost[KDIR] * mapNVars;

  #if MHD == YES
  // BX1s
  bufferSizeX3 += (data->np_tot[IDIR]+1) * data->np_tot[JDIR] * data->nghost[KDIR];
  // BX2s
  bufferSizeX3 += data->np_tot[IDIR] * (data->np_tot[JDIR]+1) * data->nghost[KDIR];
  // BX3s (not needed because reconstructed)
  bufferSizeX3 += data->np_tot[IDIR] * data->np_tot[JDIR];
  #endif  // MHD

  BufferRecvX3[faceLeft ] = IdefixArray1D<real>("BufferRecvX3Left", bufferSizeX3);
  BufferRecvX3[faceRight] = IdefixArray1D<real>("BufferRecvX3Right",bufferSizeX3);
  BufferSendX3[faceLeft ] = IdefixArray1D<real>("BufferSendX3Left", bufferSizeX3);
  BufferSendX3[faceRight] = IdefixArray1D<real>("BufferSendX3Right",bufferSizeX3);
#endif // DIMENSIONS

#ifdef MPI_PERSISTENT
  // Init persistent MPI communications
  int procSend, procRecv;

  // X1-dir exchanges
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,0,1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Send_init(BufferSendX1[faceRight].data(), bufferSizeX1, realMPI, procSend, 100,
                mygrid->CartComm, &sendRequestX1[faceRight]));

  MPI_SAFE_CALL(MPI_Recv_init(BufferRecvX1[faceLeft].data(), bufferSizeX1, realMPI, procRecv, 100,
                mygrid->CartComm, &recvRequestX1[faceLeft]));

  // Send to the left
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,0,-1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Send_init(BufferSendX1[faceLeft].data(), bufferSizeX1, realMPI, procSend, 101,
                mygrid->CartComm, &sendRequestX1[faceLeft]));

  MPI_SAFE_CALL(MPI_Recv_init(BufferRecvX1[faceRight].data(), bufferSizeX1, realMPI, procRecv, 101,
                mygrid->CartComm, &recvRequestX1[faceRight]));

  #if DIMENSIONS >= 2
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,1,1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Send_init(BufferSendX2[faceRight].data(), bufferSizeX2, realMPI, procSend, 200,
                mygrid->CartComm, &sendRequestX2[faceRight]));

  MPI_SAFE_CALL(MPI_Recv_init(BufferRecvX2[faceLeft].data(), bufferSizeX2, realMPI, procRecv, 200,
                mygrid->CartComm, &recvRequestX2[faceLeft]));

  // Send to the left
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,1,-1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Send_init(BufferSendX2[faceLeft].data(), bufferSizeX2, realMPI, procSend, 201,
                mygrid->CartComm, &sendRequestX2[faceLeft]));

  MPI_SAFE_CALL(MPI_Recv_init(BufferRecvX2[faceRight].data(), bufferSizeX2, realMPI, procRecv, 201,
                mygrid->CartComm, &recvRequestX2[faceRight]));
  #endif

  #if DIMENSIONS == 3
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,2,1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Send_init(BufferSendX3[faceRight].data(), bufferSizeX3, realMPI, procSend, 300,
                mygrid->CartComm, &sendRequestX3[faceRight]));

  MPI_SAFE_CALL(MPI_Recv_init(BufferRecvX3[faceLeft].data(), bufferSizeX3, realMPI, procRecv, 300,
                mygrid->CartComm, &recvRequestX3[faceLeft]));

  // Send to the left
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,2,-1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Send_init(BufferSendX3[faceLeft].data(), bufferSizeX3, realMPI, procSend, 301,
                mygrid->CartComm, &sendRequestX3[faceLeft]));

  MPI_SAFE_CALL(MPI_Recv_init(BufferRecvX3[faceRight].data(), bufferSizeX3, realMPI, procRecv, 301,
                mygrid->CartComm, &recvRequestX3[faceRight]));
  #endif

#endif
}

// Destructor (clean up persistent communication channels)
Mpi::~Mpi() {
  // Properly clean up the mess
#ifdef MPI_PERSISTENT
  idfx::cout << "Mpi:: Cleaning up MPI persistent communication channels" << std::endl;
  for(int i=0 ; i< 2; i++) {
    MPI_Request_free( &sendRequestX1[i]);
    MPI_Request_free( &recvRequestX1[i]);

  #if DIMENSIONS >= 2
    MPI_Request_free( &sendRequestX2[i]);
    MPI_Request_free( &recvRequestX2[i]);
  #endif

  #if DIMENSIONS == 3
    MPI_Request_free( &sendRequestX3[i]);
    MPI_Request_free( &recvRequestX3[i]);
  #endif
  }

#endif
}

void Mpi::ExchangeX1() {
  idfx::pushRegion("Mpi::ExchangeX1");

  // init MPI Timer
  idfx::mpiTimer -= timer.seconds();

  // Load  the buffers with data
  int ibeg,iend,jbeg,jend,kbeg,kend,offset;
  int nx,ny,nz;
  IdefixArray1D<real> BufferLeft=BufferSendX1[faceLeft];
  IdefixArray1D<real> BufferRight=BufferSendX1[faceRight];
  IdefixArray4D<real> Vc=data->hydro.Vc;
  IdefixArray1D<int> map = this->mapVars;
#if MHD==YES
  IdefixArray4D<real> Vs=data->hydro.Vs;
  int VsIndex;
#endif

  BoundaryType lbound = data->lbound[IDIR];
  BoundaryType rbound = data->rbound[IDIR];

  // If MPI Persistent, start receiving even before the buffers are filled
#ifdef MPI_PERSISTENT
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];

  MPI_SAFE_CALL(MPI_Startall(2, recvRequestX1));
#endif

  // Coordinates of the ghost region which needs to be transfered
  ibeg   = 0;
  iend   = data->nghost[IDIR];
  nx     = data->nghost[IDIR];  // Number of points in x
  offset = data->end[IDIR];     // Distance between beginning of left and right ghosts
  jbeg   = data->beg[JDIR];
  jend   = data->end[JDIR];
  ny     = jend - jbeg;
  kbeg   = data->beg[KDIR];
  kend   = data->end[KDIR];
  nz     = kend - kbeg;

  idefix_for("LoadBufferX1Vc",0,mapNVars,kbeg,kend,jbeg,jend,ibeg,iend,
    KOKKOS_LAMBDA (int n, int k, int j, int i) {
      BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) = Vc(map(n),k,j,i+nx);
      BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) = Vc(map(n),k,j,i+offset-nx);
    }
  );

#if MHD==YES
  // Load face-centered field in the buffer

  VsIndex = mapNVars*nx*ny*nz;

  idefix_for("LoadBufferX1BX1s",kbeg,kend,jbeg,jend,ibeg,ibeg+1,
    KOKKOS_LAMBDA (int k, int j, int i) {
      BufferLeft(i + (j-jbeg) + (k-kbeg)*ny + VsIndex ) = Vs(BX1s,k,j,i+nx);
      BufferRight(i + (j-jbeg) + (k-kbeg)*ny + VsIndex ) = Vs(BX1s,k,j,i+offset);
    }
  );

  #if DIMENSIONS >= 2
  VsIndex = mapNVars*nx*ny*nz + ny*nz;

  idefix_for("LoadBufferX1BX2s",kbeg,kend,jbeg,jend+1,ibeg,iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex ) = Vs(BX2s,k,j,i+nx);
      BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex ) = Vs(BX2s,k,j,i+offset-nx);
    }
  );

  #endif

  #if DIMENSIONS == 3
  VsIndex = mapNVars*nx*ny*nz + ny*nz + nx*(ny+1)*nz;

  idefix_for("LoadBufferX1BX3s",kbeg,kend+1,jbeg,jend,ibeg,iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(BX3s,k,j,i+nx);
      BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(BX3s,k,j,i+offset-nx);
    }
  );

  #endif
#endif

  // Wait for completion before sending out everything
  Kokkos::fence();

#ifdef MPI_PERSISTENT
  MPI_SAFE_CALL(MPI_Startall(2, sendRequestX1));
  // Wait for buffers to be received
  MPI_Waitall(2,recvRequestX1,recvStatus);

#else
  int procSend, procRecv;

  #ifdef MPI_NON_BLOCKING
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];
  MPI_Request sendRequest[2];
  MPI_Request recvRequest[2];

  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,0,1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Isend(BufferSendX1[faceRight].data(), bufferSizeX1, realMPI, procSend, 100,
                mygrid->CartComm, &sendRequest[0]));

  MPI_SAFE_CALL(MPI_Irecv(BufferRecvX1[faceLeft].data(), bufferSizeX1, realMPI, procRecv, 100,
                mygrid->CartComm, &recvRequest[0]));

  // Send to the left
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,0,-1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Isend(BufferSendX1[faceLeft].data(), bufferSizeX1, realMPI, procSend, 101,
                mygrid->CartComm, &sendRequest[1]));

  MPI_SAFE_CALL(MPI_Irecv(BufferRecvX1[faceRight].data(), bufferSizeX1, realMPI, procRecv, 101,
                mygrid->CartComm, &recvRequest[1]));

  // Wait for recv to complete (we don't care about the sends)
  MPI_Waitall(2, recvRequest, recvStatus);

  #else
  MPI_Status status;
  // Send to the right
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,0,1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX1[faceRight].data(), bufferSizeX1, realMPI, procSend, 100,
                BufferRecvX1[faceLeft].data(), bufferSizeX1, realMPI, procRecv, 100,
                mygrid->CartComm, &status));

  // Send to the left
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,0,-1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX1[faceLeft].data(), bufferSizeX1, realMPI, procSend, 101,
                BufferRecvX1[faceRight].data(), bufferSizeX1, realMPI, procRecv, 101,
                mygrid->CartComm, &status));
  #endif
#endif

  // Unpack
  BufferLeft=BufferRecvX1[faceLeft];
  BufferRight=BufferRecvX1[faceRight];

  // We fill the ghost zones
  idefix_for("StoreBufferX1Vc",0,mapNVars,kbeg,kend,jbeg,jend,ibeg,iend,
    KOKKOS_LAMBDA (int n, int k, int j, int i) {
      Vc(map(n),k,j,i) = BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz );
      Vc(map(n),k,j,i+offset) = BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz );
    }
  );

#if MHD==YES
  // Load face-centered field in the buffer

  VsIndex = mapNVars*nx*ny*nz;

  // This is not strictly required since BX1s on the neighboring domain faces
  // Should be equal by construction, since they are evolved with the same EMFs
  // HOWEVER, because the arithmetic operation are not performed in the same
  // order on the two domains, these might slowly diverge on both sides because
  // of round off error. The following loop average the two values on the faces,
  // When applicable, to eliminate this divergence between the domains
  /*
  idefix_for("StoreBufferX1BX1s",kbeg,kend,jbeg,jend,ibeg,ibeg+1,
    KOKKOS_LAMBDA (int k, int j, int i) {
      if(lbound == internal || lbound == periodic) {
        Vs(BX1s,k,j,i+nx) = HALF_F*(
                              BufferLeft(i + (j-jbeg) + (k-kbeg)*ny + VsIndex )
                            + Vs(BX1s,k,j,i+nx));
      }
      if(rbound == internal || rbound == periodic) {
        Vs(BX1s,k,j,i+offset) = HALF_F*(
                              BufferRight(i + (j-jbeg) + (k-kbeg)*ny + VsIndex )
                            + Vs(BX1s,k,j,i+offset));
      }
    }
  );
*/
  #if DIMENSIONS >= 2
  VsIndex = mapNVars*nx*ny*nz + ny*nz;

  idefix_for("StoreBufferX1BX2s",kbeg,kend,jbeg,jend+1,ibeg,iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      Vs(BX2s,k,j,i) = BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex );
      Vs(BX2s,k,j,i+offset) = BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex );
    }
  );
  #endif

  #if DIMENSIONS == 3
  VsIndex = mapNVars*nx*ny*nz + ny*nz + nx*(ny+1)*nz;

  idefix_for("StoreBufferX1BX3s",kbeg,kend+1,jbeg,jend,ibeg,iend,
    KOKKOS_LAMBDA ( int k, int j, int i) {
      Vs(BX3s,k,j,i) = BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
      Vs(BX3s,k,j,i+offset) = BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
    }
  );
  #endif
#endif

#ifdef MPI_NON_BLOCKING
  // Wait for the sends if they have not yet completed
  MPI_Waitall(2, sendRequest, sendStatus);
#endif

#ifdef MPI_PERSISTENT
  MPI_Waitall(2, sendRequestX1, sendStatus);
#endif

  // Stop MPI Timer
  idfx::mpiTimer += timer.seconds();

  idfx::popRegion();
}


void Mpi::ExchangeX2() {
  idfx::pushRegion("Mpi::ExchangeX2");

  // init MPI Timer
  idfx::mpiTimer -= timer.seconds();
  // Load  the buffers with data
  int ibeg,iend,jbeg,jend,kbeg,kend,offset;
  int nx,ny,nz;
  IdefixArray1D<real> BufferLeft=BufferSendX2[faceLeft];
  IdefixArray1D<real> BufferRight=BufferSendX2[faceRight];
  IdefixArray4D<real> Vc=data->hydro.Vc;
  IdefixArray1D<int> map = this->mapVars;
#if MHD==YES
  IdefixArray4D<real> Vs=data->hydro.Vs;
  int VsIndex;
#endif

  BoundaryType lbound = data->lbound[JDIR];
  BoundaryType rbound = data->rbound[JDIR];
// If MPI Persistent, start receiving even before the buffers are filled
#ifdef MPI_PERSISTENT
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];

  MPI_SAFE_CALL(MPI_Startall(2, recvRequestX2));
#endif

  // Coordinates of the ghost region which needs to be transfered
  ibeg   = 0;
  iend   = data->np_tot[IDIR];
  nx     = data->np_tot[IDIR];  // Number of points in x
  jbeg   = 0;
  jend   = data->nghost[JDIR];
  offset = data->end[JDIR];     // Distance between beginning of left and right ghosts
  ny     = data->nghost[JDIR];
  kbeg   = data->beg[KDIR];
  kend   = data->end[KDIR];
  nz     = kend - kbeg;

  idefix_for("LoadBufferX2Vc",0,mapNVars,kbeg,kend,jbeg,jend,ibeg,iend,
    KOKKOS_LAMBDA (int n, int k, int j, int i) {
      BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) = Vc(map(n),k,j+ny,i);
      BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) = Vc(map(n),k,j+offset-ny,i);
    }
  );

#if MHD==YES
  // Load face-centered field in the buffer

  VsIndex = mapNVars*nx*ny*nz;

  idefix_for("LoadBufferX2BX1s",kbeg,kend,jbeg,jend,ibeg,iend+1,
    KOKKOS_LAMBDA (int k, int j, int i) {
      BufferLeft(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex ) = Vs(BX1s,k,j+ny,i);
      BufferRight(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex ) = Vs(BX1s,k,j+offset-ny,i);
    }
  );

  #if DIMENSIONS >= 2
  VsIndex = mapNVars*nx*ny*nz + (nx+1)*ny*nz;

  idefix_for("LoadBufferX2BX2s",kbeg,kend,jbeg,jbeg+1,ibeg,iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx + VsIndex ) = Vs(BX2s,k,j+ny,i);
      BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx + VsIndex ) = Vs(BX2s,k,j+offset,i);
    }
  );

  #endif

  #if DIMENSIONS == 3
  VsIndex = mapNVars*nx*ny*nz + (nx+1)*ny*nz + nx*nz;

  idefix_for("LoadBufferX2BX3s",kbeg,kend+1,jbeg,jend,ibeg,iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(BX3s,k,j+ny,i);
      BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(BX3s,k,j+offset-ny,i);
    }
  );

  #endif
#endif

  // Send to the right
  Kokkos::fence();

#ifdef MPI_PERSISTENT
  MPI_SAFE_CALL(MPI_Startall(2, sendRequestX2));
  MPI_Waitall(2,recvRequestX2,recvStatus);

#else
  int procSend, procRecv;

  #ifdef MPI_NON_BLOCKING
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];
  MPI_Request sendRequest[2];
  MPI_Request recvRequest[2];

  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,1,1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Isend(BufferSendX2[faceRight].data(), bufferSizeX2, realMPI, procSend, 100,
                mygrid->CartComm, &sendRequest[0]));

  MPI_SAFE_CALL(MPI_Irecv(BufferRecvX2[faceLeft].data(), bufferSizeX2, realMPI, procRecv, 100,
                mygrid->CartComm, &recvRequest[0]));

  // Send to the left
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,1,-1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Isend(BufferSendX2[faceLeft].data(), bufferSizeX2, realMPI, procSend, 101,
                mygrid->CartComm, &sendRequest[1]));

  MPI_SAFE_CALL(MPI_Irecv(BufferRecvX2[faceRight].data(), bufferSizeX2, realMPI, procRecv, 101,
                mygrid->CartComm, &recvRequest[1]));

  // Wait for recv to complete (we don't care about the sends)
  MPI_Waitall(2, recvRequest, recvStatus);

  #else
  MPI_Status status;
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,1,1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX2[faceRight].data(), bufferSizeX2, realMPI, procSend, 200,
                BufferRecvX2[faceLeft].data(), bufferSizeX2, realMPI, procRecv, 200,
                mygrid->CartComm, &status));

  // Send to the left
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,1,-1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX2[faceLeft].data(), bufferSizeX2, realMPI, procSend, 201,
                BufferRecvX2[faceRight].data(), bufferSizeX2, realMPI, procRecv, 201,
                mygrid->CartComm, &status));
  #endif
#endif

  // Unpack
  BufferLeft=BufferRecvX2[faceLeft];
  BufferRight=BufferRecvX2[faceRight];

  // We fill the ghost zones

  idefix_for("StoreBufferX2Vc",0,mapNVars,kbeg,kend,jbeg,jend,ibeg,iend,
    KOKKOS_LAMBDA (int n, int k, int j, int i) {
      Vc(map(n),k,j,i) = BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz );
      Vc(map(n),k,j+offset,i) = BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz );
    }
  );

#if MHD==YES
  // Load face-centered field in the buffer

  VsIndex = mapNVars*nx*ny*nz;

  idefix_for("StoreBufferX2BX1s",kbeg,kend,jbeg,jend,ibeg,iend+1,
    KOKKOS_LAMBDA (int k, int j, int i) {
      Vs(BX1s,k,j,i) = BufferLeft(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex );
      Vs(BX1s,k,j+offset,i) = BufferRight(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex );
    }
  );

  #if DIMENSIONS >= 2
  VsIndex = mapNVars*nx*ny*nz + (nx+1)*ny*nz;

  idefix_for("StoreBufferX2BX2s",kbeg,kend,jbeg,jbeg+1,ibeg,iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      if(lbound == internal || lbound == periodic) {
        Vs(BX2s,k,j+ny,i) = HALF_F*(
                              BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx + VsIndex )
                            + Vs(BX2s,k,j+ny,i));
      }
      if(rbound == internal || rbound == periodic) {
        Vs(BX2s,k,j+offset,i) = HALF_F*(
                              BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx + VsIndex )
                            + Vs(BX2s,k,j+offset,i));
      }
    }
  );
  #endif

  #if DIMENSIONS == 3
  VsIndex = mapNVars*nx*ny*nz + (nx+1)*ny*nz + nx*nz;

  idefix_for("StoreBufferX2BX3s",kbeg,kend+1,jbeg,jend,ibeg,iend,
    KOKKOS_LAMBDA ( int k, int j, int i) {
      Vs(BX3s,k,j,i) = BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
      Vs(BX3s,k,j+offset,i) = BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
    }
  );
  #endif
#endif

#ifdef MPI_NON_BLOCKING
  // Wait for the sends if they have not yet completed
  MPI_Waitall(2, sendRequest, sendStatus);
#endif

#ifdef MPI_PERSISTENT
  MPI_Waitall(2, sendRequestX2, sendStatus);
#endif

  // Stop MPI Timer
  idfx::mpiTimer += timer.seconds();

  idfx::popRegion();
}


void Mpi::ExchangeX3() {
  idfx::pushRegion("Mpi::ExchangeX3");

  // init MPI Timer
  idfx::mpiTimer -= timer.seconds();

  // Load  the buffers with data
  int ibeg,iend,jbeg,jend,kbeg,kend,offset;
  int nx,ny,nz;
  IdefixArray1D<real> BufferLeft=BufferSendX3[faceLeft];
  IdefixArray1D<real> BufferRight=BufferSendX3[faceRight];
  IdefixArray4D<real> Vc=data->hydro.Vc;
  IdefixArray1D<int> map = this->mapVars;
#if MHD==YES
  IdefixArray4D<real> Vs=data->hydro.Vs;
  int VsIndex;
#endif

  BoundaryType lbound = data->lbound[KDIR];
  BoundaryType rbound = data->rbound[KDIR];

  // If MPI Persistent, start receiving even before the buffers are filled
#ifdef MPI_PERSISTENT
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];

  MPI_SAFE_CALL(MPI_Startall(2, recvRequestX3));
#endif

  // Coordinates of the ghost region which needs to be transfered
  ibeg   = 0;
  iend   = data->np_tot[IDIR];
  nx     = data->np_tot[IDIR];  // Number of points in x
  jbeg   = 0;
  jend   = data->np_tot[JDIR];
  ny     = data->np_tot[JDIR];
  kbeg   = 0;
  kend   = data->nghost[KDIR];
  offset = data->end[KDIR];     // Distance between beginning of left and right ghosts
  nz     = data->nghost[KDIR];

  idefix_for("LoadBufferX3Vc",0,mapNVars,kbeg,kend,jbeg,jend,ibeg,iend,
      KOKKOS_LAMBDA (int n, int k, int j, int i) {
        int nt = map(n);
        BufferLeft((i-ibeg) + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) = Vc(nt,k+nz,j,i);
        BufferRight((i-ibeg) + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) = Vc(nt,k+offset-nz,j,i);
      }
  );

#if MHD==YES
  // Load face-centered field in the buffer

  VsIndex = mapNVars*nx*ny*nz;

  idefix_for("LoadBufferX3BX1s",kbeg,kend,jbeg,jend,ibeg,iend+1,
    KOKKOS_LAMBDA (int k, int j, int i) {
      BufferLeft(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex ) = Vs(BX1s,k+nz,j,i);
      BufferRight(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex ) = Vs(BX1s,k+offset-nz,j,i);
    }
  );

  #if DIMENSIONS >= 2
  VsIndex = mapNVars*nx*ny*nz + (nx+1)*ny*nz;

  idefix_for("LoadBufferX3BX2s",kbeg,kend,jbeg,jend+1,ibeg,iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex ) = Vs(BX2s,k+nz,j,i);
      BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex ) = Vs(BX2s,k+offset-nz,j,i);
    }
  );
  #endif

  #if DIMENSIONS == 3
  VsIndex = mapNVars*nx*ny*nz + (nx+1)*ny*nz + nx*(ny+1)*nz;

  idefix_for("LoadBufferX3BX3s",kbeg,kbeg+1,jbeg,jend,ibeg,iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(BX3s,k+nz,j,i);
      BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(BX3s,k+offset,j,i);
    }
  );
  #endif
#endif

  // Send to the right
  Kokkos::fence();

#ifdef MPI_PERSISTENT
  MPI_SAFE_CALL(MPI_Startall(2, sendRequestX3));
  MPI_Waitall(2,recvRequestX3,recvStatus);

#else
  int procSend, procRecv;

  #ifdef MPI_NON_BLOCKING
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];
  MPI_Request sendRequest[2];
  MPI_Request recvRequest[2];

  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,2,1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Isend(BufferSendX3[faceRight].data(), bufferSizeX3, realMPI, procSend, 100,
                mygrid->CartComm, &sendRequest[0]));

  MPI_SAFE_CALL(MPI_Irecv(BufferRecvX3[faceLeft].data(), bufferSizeX3, realMPI, procRecv, 100,
                mygrid->CartComm, &recvRequest[0]));

  // Send to the left
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,2,-1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Isend(BufferSendX3[faceLeft].data(), bufferSizeX3, realMPI, procSend, 101,
                mygrid->CartComm, &sendRequest[1]));

  MPI_SAFE_CALL(MPI_Irecv(BufferRecvX3[faceRight].data(), bufferSizeX3, realMPI, procRecv, 101,
                mygrid->CartComm, &recvRequest[1]));

  // Wait for recv to complete (we don't care about the sends)
  MPI_Waitall(2, recvRequest, recvStatus);

  #else
  MPI_Status status;
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,2,1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX3[faceRight].data(), bufferSizeX3, realMPI, procSend, 300,
                BufferRecvX3[faceLeft].data(), bufferSizeX3, realMPI, procRecv, 300,
                mygrid->CartComm, &status));

  // Send to the left
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,2,-1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX3[faceLeft].data(), bufferSizeX3, realMPI, procSend, 301,
                BufferRecvX3[faceRight].data(), bufferSizeX3, realMPI, procRecv, 301,
                mygrid->CartComm, &status));
  #endif
#endif

  // Unpack
  BufferLeft=BufferRecvX3[faceLeft];
  BufferRight=BufferRecvX3[faceRight];

  // We fill the ghost zones
  idefix_for("StoreBufferX3Vc",0,mapNVars,kbeg,kend,jbeg,jend,ibeg,iend,
    KOKKOS_LAMBDA (int n, int k, int j, int i) {
      Vc(map(n),k,j,i) = BufferLeft((i-ibeg) + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz );
      Vc(map(n),k+offset,j,i) = BufferRight((i-ibeg) + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz );
    }
  );

#if MHD==YES
  // Load face-centered field in the buffer

  VsIndex = mapNVars*nx*ny*nz;

  idefix_for("StoreBufferX3BX1s",kbeg,kend,jbeg,jend,ibeg,iend+1,
    KOKKOS_LAMBDA (int k, int j, int i) {
      Vs(BX1s,k,j,i) = BufferLeft(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex );
      Vs(BX1s,k+offset,j,i) = BufferRight(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex );
    }
  );

  #if DIMENSIONS >= 2
  VsIndex = mapNVars*nx*ny*nz + (nx+1)*ny*nz;

  idefix_for("StoreBufferX3BX2s",kbeg,kend,jbeg,jend+1,ibeg,iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      Vs(BX2s,k,j,i) = BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex );
      Vs(BX2s,k+offset,j,i) = BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex );
    }
  );
  #endif

  #if DIMENSIONS == 3
  VsIndex = mapNVars*nx*ny*nz + (nx+1)*ny*nz + nx*(ny+1)*nz;

  idefix_for("StoreBufferX3BX3s",kbeg,kbeg+1,jbeg,jend,ibeg,iend,
    KOKKOS_LAMBDA ( int k, int j, int i) {
      if(lbound == internal || lbound == periodic) {
        Vs(BX3s,k+nz,j,i) = HALF_F*(
                              BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex )
                            + Vs(BX3s,k+nz,j,i));
      }
      if(rbound == internal || rbound == periodic) {
        Vs(BX3s,k+offset,j,i) = HALF_F*(
                                    BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex )
                                  + Vs(BX3s,k+offset,j,i));
      }
    }
  );
  #endif

#endif

#ifdef MPI_NON_BLOCKING
  // Wait for the sends if they have not yet completed
  MPI_Waitall(2, sendRequest, sendStatus);
#endif

#ifdef MPI_PERSISTENT
  MPI_Waitall(2, sendRequestX3, sendStatus);
#endif

  // Stop MPI Timer
  idfx::mpiTimer += timer.seconds();

  idfx::popRegion();
}
