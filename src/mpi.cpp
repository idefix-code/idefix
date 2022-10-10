// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************


#include "mpi.hpp"
#include <signal.h>
#include <string>
#include "idefix.hpp"
#include "dataBlock.hpp"


#if defined(OPEN_MPI) && OPEN_MPI
#include "mpi-ext.h"                // Needed for CUDA-aware check */
#endif


//#define MPI_NON_BLOCKING
#define MPI_PERSISTENT

// init the number of instances
int Mpi::nInstances = 0;

// MPI Routines exchange
void Mpi::ExchangeAll() {
  IDEFIX_ERROR("Not Implemented");
}

///
/// Initialise an instance of the MPI class.
/// @param grid: pointer to the grid object (needed to get the MPI neighbours)
/// @param inputMap: 1st indices of inputVc which are to be exchanged (i.e, the list of variables)
/// @param nghost: size of the ghost region in each direction
/// @param nint: size of the internal region in each direction
/// @param inputHaveVs: whether the instance should also treat face-centered variable
///                     (optional, default false)
///

void Mpi::Init(Grid *grid, std::vector<int> inputMap,
               int nghost[3], int nint[3],
               bool inputHaveVs) {
  this->mygrid = grid;

  // increase the number of instances
  nInstances++;
  thisInstance=nInstances;

  // Transfer the vector of indices as an IdefixArray on the target

  // Allocate mapVars on target and copy it from the input argument list
  this->mapVars = idfx::ConvertVectorToIdefixArray(inputMap);
  this->mapNVars = inputMap.size();
  this->haveVs = inputHaveVs;

  // Compute indices of arrays we will be working with
  for(int dir = 0 ; dir < 3 ; dir++) {
    this->nghost[dir] = nghost[dir];
    this->nint[dir] = nint[dir];
    this->ntot[dir] = nint[dir]+2*nghost[dir];
    this->beg[dir] = nghost[dir];
    this->end[dir] = nghost[dir]+nint[dir];
  }

  /////////////////////////////////////////////////////////////////////////////
  // Init exchange datasets
  bufferSizeX1 = 0;
  bufferSizeX2 = 0;
  bufferSizeX3 = 0;

  // Number of cells in X1 boundary condition:
  bufferSizeX1 = nghost[IDIR] * nint[JDIR] * nint[KDIR] * mapNVars;

  if(haveVs) {
    #if DIMENSIONS>=2
    bufferSizeX1 += nghost[IDIR] * (nint[JDIR]+1) * nint[KDIR];
    #endif

    #if DIMENSIONS==3
    bufferSizeX1 += nghost[IDIR] * nint[JDIR] * (nint[KDIR]+1);
    #endif  // DIMENSIONS
  }


  BufferRecvX1[faceLeft ] = IdefixArray1D<real>("BufferRecvX1Left", bufferSizeX1);
  BufferRecvX1[faceRight] = IdefixArray1D<real>("BufferRecvX1Right",bufferSizeX1);
  BufferSendX1[faceLeft ] = IdefixArray1D<real>("BufferSendX1Left", bufferSizeX1);
  BufferSendX1[faceRight] = IdefixArray1D<real>("BufferSendX1Right",bufferSizeX1);

  // Number of cells in X2 boundary condition (only required when problem >2D):
#if DIMENSIONS >= 2
  bufferSizeX2 = ntot[IDIR] * nghost[JDIR] * nint[KDIR] * mapNVars;
  if(haveVs) {
    // IDIR
    bufferSizeX2 += (ntot[IDIR]+1) * nghost[JDIR] * nint[KDIR];
    #if DIMENSIONS==3
    bufferSizeX2 += ntot[IDIR] * nghost[JDIR] * (nint[KDIR]+1);
    #endif  // DIMENSIONS
  }

  BufferRecvX2[faceLeft ] = IdefixArray1D<real>("BufferRecvX2Left", bufferSizeX2);
  BufferRecvX2[faceRight] = IdefixArray1D<real>("BufferRecvX2Right",bufferSizeX2);
  BufferSendX2[faceLeft ] = IdefixArray1D<real>("BufferSendX2Left", bufferSizeX2);
  BufferSendX2[faceRight] = IdefixArray1D<real>("BufferSendX2Right",bufferSizeX2);

#endif
// Number of cells in X3 boundary condition (only required when problem is 3D):
#if DIMENSIONS ==3
  bufferSizeX3 = ntot[IDIR] * ntot[JDIR] * nghost[KDIR] * mapNVars;

  if(haveVs) {
    // IDIR
    bufferSizeX3 += (ntot[IDIR]+1) * ntot[JDIR] * nghost[KDIR];
    // JDIR
    bufferSizeX3 += ntot[IDIR] * (ntot[JDIR]+1) * nghost[KDIR];
  }

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

  MPI_SAFE_CALL(MPI_Send_init(BufferSendX1[faceRight].data(), bufferSizeX1, realMPI, procSend,
                thisInstance*1000, mygrid->CartComm, &sendRequestX1[faceRight]));

  MPI_SAFE_CALL(MPI_Recv_init(BufferRecvX1[faceLeft].data(), bufferSizeX1, realMPI, procRecv,
                thisInstance*1000, mygrid->CartComm, &recvRequestX1[faceLeft]));

  // Send to the left
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,0,-1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Send_init(BufferSendX1[faceLeft].data(), bufferSizeX1, realMPI, procSend,
                thisInstance*1000+1,mygrid->CartComm, &sendRequestX1[faceLeft]));

  MPI_SAFE_CALL(MPI_Recv_init(BufferRecvX1[faceRight].data(), bufferSizeX1, realMPI, procRecv,
                thisInstance*1000+1,mygrid->CartComm, &recvRequestX1[faceRight]));

  #if DIMENSIONS >= 2
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,1,1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Send_init(BufferSendX2[faceRight].data(), bufferSizeX2, realMPI, procSend,
                thisInstance*1000+10, mygrid->CartComm, &sendRequestX2[faceRight]));

  MPI_SAFE_CALL(MPI_Recv_init(BufferRecvX2[faceLeft].data(), bufferSizeX2, realMPI, procRecv,
                thisInstance*1000+10, mygrid->CartComm, &recvRequestX2[faceLeft]));

  // Send to the left
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,1,-1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Send_init(BufferSendX2[faceLeft].data(), bufferSizeX2, realMPI, procSend,
                thisInstance*1000+11, mygrid->CartComm, &sendRequestX2[faceLeft]));

  MPI_SAFE_CALL(MPI_Recv_init(BufferRecvX2[faceRight].data(), bufferSizeX2, realMPI, procRecv,
                thisInstance*1000+11, mygrid->CartComm, &recvRequestX2[faceRight]));
  #endif

  #if DIMENSIONS == 3
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,2,1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Send_init(BufferSendX3[faceRight].data(), bufferSizeX3, realMPI, procSend,
                thisInstance*1000+20, mygrid->CartComm, &sendRequestX3[faceRight]));

  MPI_SAFE_CALL(MPI_Recv_init(BufferRecvX3[faceLeft].data(), bufferSizeX3, realMPI, procRecv,
                thisInstance*1000+20, mygrid->CartComm, &recvRequestX3[faceLeft]));

  // Send to the left
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,2,-1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Send_init(BufferSendX3[faceLeft].data(), bufferSizeX3, realMPI, procSend,
                thisInstance*1000+21, mygrid->CartComm, &sendRequestX3[faceLeft]));

  MPI_SAFE_CALL(MPI_Recv_init(BufferRecvX3[faceRight].data(), bufferSizeX3, realMPI, procRecv,
                thisInstance*1000+21, mygrid->CartComm, &recvRequestX3[faceRight]));
  #endif

#endif // MPI_Persistent

  // say this instance is initialized.
  isInitialized = true;
}

// Destructor (clean up persistent communication channels)
Mpi::~Mpi() {
  if(isInitialized) {
    // Properly clean up the mess
    #ifdef MPI_PERSISTENT
      idfx::cout << "Mpi(" << thisInstance
                << "): Cleaning up MPI persistent communication channels" << std::endl;
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
    if(thisInstance==1) {
      idfx::cout << "Mpi(" << thisInstance << "): measured throughput is "
                << bytesSentOrReceived/myTimer/1024.0/1024.0 << " MB/s" << std::endl;
      idfx::cout << "Mpi(" << thisInstance << "): message sizes were " << std::endl;
      idfx::cout << "        X1: " << bufferSizeX1*sizeof(real)/1024.0/1024.0 << " MB" << std::endl;
      idfx::cout << "        X2: " << bufferSizeX2*sizeof(real)/1024.0/1024.0 << " MB" << std::endl;
      idfx::cout << "        X3: " << bufferSizeX3*sizeof(real)/1024.0/1024.0 << " MB" << std::endl;
    }
    isInitialized = false;
  }
}

void Mpi::ExchangeX1(IdefixArray4D<real> Vc, IdefixArray4D<real> Vs) {
  idfx::pushRegion("Mpi::ExchangeX1");

  // Load  the buffers with data
  int ibeg,iend,jbeg,jend,kbeg,kend,offset;
  int nx,ny,nz;
  IdefixArray1D<real> BufferLeft=BufferSendX1[faceLeft];
  IdefixArray1D<real> BufferRight=BufferSendX1[faceRight];
  IdefixArray1D<int> map = this->mapVars;

  // If MPI Persistent, start receiving even before the buffers are filled
  myTimer -= MPI_Wtime();
#ifdef MPI_PERSISTENT
  double tStart = MPI_Wtime();
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];

  MPI_SAFE_CALL(MPI_Startall(2, recvRequestX1));
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;
#endif
  myTimer += MPI_Wtime();

  // Coordinates of the ghost region which needs to be transfered
  ibeg   = 0;
  iend   = nghost[IDIR];
  nx     = nghost[IDIR];  // Number of points in x
  offset = end[IDIR];     // Distance between beginning of left and right ghosts
  jbeg   = beg[JDIR];
  jend   = end[JDIR];
  ny     = jend - jbeg;
  kbeg   = beg[KDIR];
  kend   = end[KDIR];
  nz     = kend - kbeg;

  idefix_for("LoadBufferX1Vc",0,mapNVars,kbeg,kend,jbeg,jend,ibeg,iend,
    KOKKOS_LAMBDA (int n, int k, int j, int i) {
      BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) = Vc(map(n),k,j,i+nx);
      BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) = Vc(map(n),k,j,i+offset-nx);
    }
  );

  // Load face-centered field in the buffer
  if(haveVs) {
    #if DIMENSIONS >= 2
    int VsIndex = mapNVars*nx*ny*nz;

    idefix_for("LoadBufferX1JDIR",kbeg,kend,jbeg,jend+1,ibeg,iend,
      KOKKOS_LAMBDA (int k, int j, int i) {
        BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex ) = Vs(JDIR,k,j,i+nx);
        BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex ) = Vs(JDIR,k,j,i+offset-nx);
      }
    );

    #endif

    #if DIMENSIONS == 3
    VsIndex = mapNVars*nx*ny*nz + nx*(ny+1)*nz;

    idefix_for("LoadBufferX1KDIR",kbeg,kend+1,jbeg,jend,ibeg,iend,
      KOKKOS_LAMBDA (int k, int j, int i) {
        BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(KDIR,k,j,i+nx);
        BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(KDIR,k,j,i+offset-nx);
      }
    );

    #endif
  }

  // Wait for completion before sending out everything
  Kokkos::fence();
  myTimer -= MPI_Wtime();
#ifdef MPI_PERSISTENT
  tStart = MPI_Wtime();
  MPI_SAFE_CALL(MPI_Startall(2, sendRequestX1));
  // Wait for buffers to be received
  MPI_Waitall(2,recvRequestX1,recvStatus);
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;

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

  tStart = MPI_Wtime();
  MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX1[faceLeft].data(), bufferSizeX1, realMPI, procSend, 101,
                BufferRecvX1[faceRight].data(), bufferSizeX1, realMPI, procRecv, 101,
                mygrid->CartComm, &status));
  #endif
#endif
  myTimer += MPI_Wtime();
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

  if(haveVs) {
    #if DIMENSIONS >= 2
    int VsIndex = mapNVars*nx*ny*nz;

    idefix_for("StoreBufferX1JDIR",kbeg,kend,jbeg,jend+1,ibeg,iend,
      KOKKOS_LAMBDA (int k, int j, int i) {
        Vs(JDIR,k,j,i) = BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex );
        Vs(JDIR,k,j,i+offset) = BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex );
      }
    );
    #endif

    #if DIMENSIONS == 3
    VsIndex = mapNVars*nx*ny*nz + nx*(ny+1)*nz;

    idefix_for("StoreBufferX1KDIR",kbeg,kend+1,jbeg,jend,ibeg,iend,
      KOKKOS_LAMBDA ( int k, int j, int i) {
        Vs(KDIR,k,j,i) = BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
        Vs(KDIR,k,j,i+offset) = BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
      }
    );
    #endif
  }

myTimer -= MPI_Wtime();
#ifdef MPI_NON_BLOCKING
  // Wait for the sends if they have not yet completed
  MPI_Waitall(2, sendRequest, sendStatus);
#endif

#ifdef MPI_PERSISTENT
  MPI_Waitall(2, sendRequestX1, sendStatus);
#endif
  myTimer += MPI_Wtime();
  bytesSentOrReceived += 4*bufferSizeX1*sizeof(real);

  idfx::popRegion();
}


void Mpi::ExchangeX2(IdefixArray4D<real> Vc, IdefixArray4D<real> Vs) {
  idfx::pushRegion("Mpi::ExchangeX2");

  // Load  the buffers with data
  int ibeg,iend,jbeg,jend,kbeg,kend,offset;
  int nx,ny,nz;
  IdefixArray1D<real> BufferLeft=BufferSendX2[faceLeft];
  IdefixArray1D<real> BufferRight=BufferSendX2[faceRight];
  IdefixArray1D<int> map = this->mapVars;

// If MPI Persistent, start receiving even before the buffers are filled
  myTimer -= MPI_Wtime();
#ifdef MPI_PERSISTENT
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];

  double tStart = MPI_Wtime();
  MPI_SAFE_CALL(MPI_Startall(2, recvRequestX2));
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;
#endif
  myTimer += MPI_Wtime();

  // Coordinates of the ghost region which needs to be transfered
  ibeg   = 0;
  iend   = ntot[IDIR];
  nx     = ntot[IDIR];  // Number of points in x
  jbeg   = 0;
  jend   = nghost[JDIR];
  offset = end[JDIR];     // Distance between beginning of left and right ghosts
  ny     = nghost[JDIR];
  kbeg   = beg[KDIR];
  kend   = end[KDIR];
  nz     = kend - kbeg;

  idefix_for("LoadBufferX2Vc",0,mapNVars,kbeg,kend,jbeg,jend,ibeg,iend,
    KOKKOS_LAMBDA (int n, int k, int j, int i) {
      BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) = Vc(map(n),k,j+ny,i);
      BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) = Vc(map(n),k,j+offset-ny,i);
    }
  );

  // Load face-centered field in the buffer
  if(haveVs) {
    int VsIndex = mapNVars*nx*ny*nz;

    idefix_for("LoadBufferX2IDIR",kbeg,kend,jbeg,jend,ibeg,iend+1,
      KOKKOS_LAMBDA (int k, int j, int i) {
        BufferLeft(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex ) = Vs(IDIR,k,j+ny,i);
        BufferRight(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex ) = Vs(IDIR,k,j+offset-ny,i);
      }
    );


    #if DIMENSIONS == 3
    VsIndex = mapNVars*nx*ny*nz + (nx+1)*ny*nz;

    idefix_for("LoadBufferX2KDIR",kbeg,kend+1,jbeg,jend,ibeg,iend,
      KOKKOS_LAMBDA (int k, int j, int i) {
        BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(KDIR,k,j+ny,i);
        BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(KDIR,k,j+offset-ny,i);
      }
    );

    #endif
  }

  // Send to the right
  Kokkos::fence();

  myTimer -= MPI_Wtime();
#ifdef MPI_PERSISTENT
  tStart = MPI_Wtime();
  MPI_SAFE_CALL(MPI_Startall(2, sendRequestX2));
  MPI_Waitall(2,recvRequestX2,recvStatus);
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;

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

  tStart = MPI_Wtime():
  MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX2[faceLeft].data(), bufferSizeX2, realMPI, procSend, 201,
                BufferRecvX2[faceRight].data(), bufferSizeX2, realMPI, procRecv, 201,
                mygrid->CartComm, &status));
  #endif
#endif
  myTimer += MPI_Wtime();
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

  // Load face-centered field in the buffer
  if(haveVs) {
    int VsIndex = mapNVars*nx*ny*nz;

    idefix_for("StoreBufferX2IDIR",kbeg,kend,jbeg,jend,ibeg,iend+1,
      KOKKOS_LAMBDA (int k, int j, int i) {
        Vs(IDIR,k,j,i) = BufferLeft(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex );
        Vs(IDIR,k,j+offset,i) = BufferRight(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex );
      }
    );

    #if DIMENSIONS == 3
    VsIndex = mapNVars*nx*ny*nz + (nx+1)*ny*nz;

    idefix_for("StoreBufferX2KDIR",kbeg,kend+1,jbeg,jend,ibeg,iend,
      KOKKOS_LAMBDA ( int k, int j, int i) {
        Vs(KDIR,k,j,i) = BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
        Vs(KDIR,k,j+offset,i) = BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
      }
    );
    #endif
  }

  myTimer -= MPI_Wtime();
#ifdef MPI_NON_BLOCKING
  // Wait for the sends if they have not yet completed
  MPI_Waitall(2, sendRequest, sendStatus);
#endif

#ifdef MPI_PERSISTENT
  MPI_Waitall(2, sendRequestX2, sendStatus);
#endif
  myTimer += MPI_Wtime();
  bytesSentOrReceived += 4*bufferSizeX2*sizeof(real);

  idfx::popRegion();
}


void Mpi::ExchangeX3(IdefixArray4D<real> Vc, IdefixArray4D<real> Vs) {
  idfx::pushRegion("Mpi::ExchangeX3");


  // Load  the buffers with data
  int ibeg,iend,jbeg,jend,kbeg,kend,offset;
  int nx,ny,nz;
  IdefixArray1D<real> BufferLeft=BufferSendX3[faceLeft];
  IdefixArray1D<real> BufferRight=BufferSendX3[faceRight];
  IdefixArray1D<int> map = this->mapVars;

  // If MPI Persistent, start receiving even before the buffers are filled
  myTimer -= MPI_Wtime();
#ifdef MPI_PERSISTENT
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];

  double tStart = MPI_Wtime();
  MPI_SAFE_CALL(MPI_Startall(2, recvRequestX3));
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;
#endif
  myTimer += MPI_Wtime();
  // Coordinates of the ghost region which needs to be transfered
  ibeg   = 0;
  iend   = ntot[IDIR];
  nx     = ntot[IDIR];  // Number of points in x
  jbeg   = 0;
  jend   = ntot[JDIR];
  ny     = ntot[JDIR];
  kbeg   = 0;
  kend   = nghost[KDIR];
  offset = end[KDIR];     // Distance between beginning of left and right ghosts
  nz     = nghost[KDIR];

  idefix_for("LoadBufferX3Vc",0,mapNVars,kbeg,kend,jbeg,jend,ibeg,iend,
      KOKKOS_LAMBDA (int n, int k, int j, int i) {
        int nt = map(n);
        BufferLeft((i-ibeg) + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) = Vc(nt,k+nz,j,i);
        BufferRight((i-ibeg) + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) = Vc(nt,k+offset-nz,j,i);
      }
  );

  // Load face-centered field in the buffer
  if(haveVs) {
    int VsIndex = mapNVars*nx*ny*nz;

    idefix_for("LoadBufferX3IDIR",kbeg,kend,jbeg,jend,ibeg,iend+1,
      KOKKOS_LAMBDA (int k, int j, int i) {
        BufferLeft(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex ) = Vs(IDIR,k+nz,j,i);
        BufferRight(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex ) = Vs(IDIR,k+offset-nz,j,i);
      }
    );

    #if DIMENSIONS >= 2
    VsIndex = mapNVars*nx*ny*nz + (nx+1)*ny*nz;

    idefix_for("LoadBufferX3JDIR",kbeg,kend,jbeg,jend+1,ibeg,iend,
      KOKKOS_LAMBDA (int k, int j, int i) {
        BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex ) = Vs(JDIR,k+nz,j,i);
        BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex ) = Vs(JDIR,k+offset-nz,j,i);
      }
    );
    #endif
  }

  // Send to the right
  Kokkos::fence();

  myTimer -= MPI_Wtime();
#ifdef MPI_PERSISTENT
  tStart = MPI_Wtime();
  MPI_SAFE_CALL(MPI_Startall(2, sendRequestX3));
  MPI_Waitall(2,recvRequestX3,recvStatus);
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;

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

  tStart = MPI_Wtime();
  MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX3[faceLeft].data(), bufferSizeX3, realMPI, procSend, 301,
                BufferRecvX3[faceRight].data(), bufferSizeX3, realMPI, procRecv, 301,
                mygrid->CartComm, &status));
  #endif
#endif
  myTimer += MPI_Wtime();
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

  // Load face-centered field in the buffer
  if(haveVs) {
    int VsIndex = mapNVars*nx*ny*nz;

    idefix_for("StoreBufferX3IDIR",kbeg,kend,jbeg,jend,ibeg,iend+1,
      KOKKOS_LAMBDA (int k, int j, int i) {
        Vs(IDIR,k,j,i) = BufferLeft(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex );
        Vs(IDIR,k+offset,j,i) = BufferRight(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex );
      }
    );

    #if DIMENSIONS >= 2
    VsIndex = mapNVars*nx*ny*nz + (nx+1)*ny*nz;

    idefix_for("StoreBufferX3JDIR",kbeg,kend,jbeg,jend+1,ibeg,iend,
      KOKKOS_LAMBDA (int k, int j, int i) {
        Vs(JDIR,k,j,i) = BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex );
        Vs(JDIR,k+offset,j,i) = BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex );
      }
    );
    #endif
  }

  myTimer -= MPI_Wtime();
#ifdef MPI_NON_BLOCKING
  // Wait for the sends if they have not yet completed
  MPI_Waitall(2, sendRequest, sendStatus);
#endif

#ifdef MPI_PERSISTENT
  MPI_Waitall(2, sendRequestX3, sendStatus);
#endif
  myTimer += MPI_Wtime();
  bytesSentOrReceived += 4*bufferSizeX2*sizeof(real);

  idfx::popRegion();
}



void Mpi::CheckConfig() {
  idfx::pushRegion("Mpi::CheckConfig");
  // compile time check
  #ifdef KOKKOS_ENABLE_CUDA
    #if defined(MPIX_CUDA_AWARE_SUPPORT) && !MPIX_CUDA_AWARE_SUPPORT
      #error Your MPI library is not CUDA Aware (check Idefix requirements).
    #endif
  #endif /* MPIX_CUDA_AWARE_SUPPORT */

  // Run-time check that we can do a reduce on device arrays
  IdefixArray1D<int64_t> src("MPIChecksrc",1);
  IdefixArray1D<int64_t>::HostMirror srcHost = Kokkos::create_mirror_view(src);

  srcHost(0) = idfx::prank;
  Kokkos::deep_copy(src, srcHost);

  IdefixArray1D<int64_t> dst("MPICheckdst",1);
  IdefixArray1D<int64_t>::HostMirror dstHost = Kokkos::create_mirror_view(dst);

  MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

  // Capture segfaults
  struct sigaction newHandler;
  struct sigaction oldHandler;
  memset(&newHandler, 0, sizeof(newHandler));
  newHandler.sa_flags = SA_SIGINFO;
  newHandler.sa_sigaction = Mpi::SigErrorHandler;
  sigaction(SIGSEGV, &newHandler, &oldHandler);

  try {
    int ierr = MPI_Allreduce(src.data(), dst.data(), 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    if(ierr != 0) {
      char MPImsg[MPI_MAX_ERROR_STRING];
      int MPImsgLen;
      MPI_Error_string(ierr, MPImsg, &MPImsgLen);
      throw std::runtime_error(std::string(MPImsg, MPImsgLen));
    }
  } catch(std::exception &e) {
    std::stringstream errmsg;
    errmsg << "Your MPI library is unable to perform reductions on Idefix arrays.";
    errmsg << std::endl;
    #ifdef KOKKOS_ENABLE_CUDA
      errmsg << "Check that your MPI library is CUDA aware." << std::endl;
    #elif KOKKOS_ENABLE_HIP
      errmsg << "Check that your MPI library is RocM aware." << std::endl;
    #else
      errmsg << "Check your MPI library configuration." << std::endl;
    #endif
    errmsg << "Error: " << e.what() << std::endl;
    IDEFIX_ERROR(errmsg);
  }
  // Restore old handlers
  sigaction(SIGSEGV, &oldHandler, NULL );
  MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);

    // Check that we have the proper end result
  Kokkos::deep_copy(dstHost, dst);
  int64_t size = static_cast<int64_t>(idfx::psize);
  if(dstHost(0) != size*(size-1)/2) {
    std::stringstream errmsg;
    errmsg << "Your MPI library managed to perform reduction on Idefix Arrays, but the result ";
    errmsg << "is incorrect. " << std::endl;
    errmsg << "Check your MPI library configuration." << std::endl;
    IDEFIX_ERROR(errmsg);
  }
  idfx::popRegion();
}

void Mpi::SigErrorHandler(int nSignum, siginfo_t* si, void* vcontext) {
  std::stringstream errmsg;
  errmsg << "A segmentation fault was triggered while attempting to test your MPI library.";
  errmsg << std::endl;
  errmsg << "Your MPI library is unable to perform reductions on Idefix arrays.";
  errmsg << std::endl;
  #ifdef KOKKOS_ENABLE_CUDA
    errmsg << "Check that your MPI library is CUDA aware." << std::endl;
  #elif KOKKOS_ENABLE_HIP
    errmsg << "Check that your MPI library is RocM aware." << std::endl;
  #else
    errmsg << "Check your MPI library configuration." << std::endl;
  #endif
  IDEFIX_ERROR(errmsg);
}
