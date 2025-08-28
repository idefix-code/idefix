// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************


#include "mpi.hpp"
#include <signal.h>
#include <string>
#include <chrono>   // NOLINT [build/c++11]
#include <thread>  // NOLINT [build/c++11]
#include <utility>
#include <vector>
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
  idfx::pushRegion("Mpi::Init");
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
    bufferSizeX1 += nghost[IDIR] * nint[JDIR] * nint[KDIR];
    #if DIMENSIONS>=2
    bufferSizeX1 += nghost[IDIR] * (nint[JDIR]+1) * nint[KDIR];
    #endif

    #if DIMENSIONS==3
    bufferSizeX1 += nghost[IDIR] * nint[JDIR] * (nint[KDIR]+1);
    #endif  // DIMENSIONS
  }


  BufferRecvX1[faceLeft ] = Buffer(bufferSizeX1);
  BufferRecvX1[faceRight] = Buffer(bufferSizeX1);
  BufferSendX1[faceLeft ] = Buffer(bufferSizeX1);
  BufferSendX1[faceRight] = Buffer(bufferSizeX1);

  // Number of cells in X2 boundary condition (only required when problem >2D):
#if DIMENSIONS >= 2
  bufferSizeX2 = ntot[IDIR] * nghost[JDIR] * nint[KDIR] * mapNVars;
  if(haveVs) {
    // IDIR
    bufferSizeX2 += (ntot[IDIR]+1) * nghost[JDIR] * nint[KDIR];
    #if DIMENSIONS>=2
    bufferSizeX2 += ntot[IDIR] * nghost[JDIR] * nint[KDIR];
    #endif
    #if DIMENSIONS==3
    bufferSizeX2 += ntot[IDIR] * nghost[JDIR] * (nint[KDIR]+1);
    #endif  // DIMENSIONS
  }

  BufferRecvX2[faceLeft ] = Buffer(bufferSizeX2);
  BufferRecvX2[faceRight] = Buffer(bufferSizeX2);
  BufferSendX2[faceLeft ] = Buffer(bufferSizeX2);
  BufferSendX2[faceRight] = Buffer(bufferSizeX2);

#endif
// Number of cells in X3 boundary condition (only required when problem is 3D):
#if DIMENSIONS ==3
  bufferSizeX3 = ntot[IDIR] * ntot[JDIR] * nghost[KDIR] * mapNVars;

  if(haveVs) {
    // IDIR
    bufferSizeX3 += (ntot[IDIR]+1) * ntot[JDIR] * nghost[KDIR];
    // JDIR
    bufferSizeX3 += ntot[IDIR] * (ntot[JDIR]+1) * nghost[KDIR];
    // KDIR
    bufferSizeX3 += ntot[IDIR] * ntot[JDIR] * nghost[KDIR];
  }

  BufferRecvX3[faceLeft ] = Buffer(bufferSizeX3);
  BufferRecvX3[faceRight] = Buffer(bufferSizeX3);
  BufferSendX3[faceLeft ] = Buffer(bufferSizeX3);
  BufferSendX3[faceRight] = Buffer(bufferSizeX3);
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

  idfx::popRegion();
}

// Destructor (clean up persistent communication channels)
Mpi::~Mpi() {
  idfx::pushRegion("Mpi::~Mpi");
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
  idfx::popRegion();
}

void Mpi::ExchangeX1(IdefixArray4D<real> Vc, IdefixArray4D<real> Vs) {
  idfx::pushRegion("Mpi::ExchangeX1");

  // Load  the buffers with data
  int ibeg,iend,jbeg,jend,kbeg,kend,offset,nx;
  Buffer BufferLeft = BufferSendX1[faceLeft];
  Buffer BufferRight = BufferSendX1[faceRight];
  IdefixArray1D<int> map = this->mapVars;

  // If MPI Persistent, start receiving even before the buffers are filled
  myTimer -= MPI_Wtime();
  double tStart = MPI_Wtime();
#ifdef MPI_PERSISTENT
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

  kbeg   = beg[KDIR];
  kend   = end[KDIR];


  BufferLeft.ResetPointer();
  BufferRight.ResetPointer();

  BufferLeft.Pack(Vc, map, std::make_pair(ibeg+nx, iend+nx),
                           std::make_pair(jbeg   , jend),
                           std::make_pair(kbeg   , kend));

  BufferRight.Pack(Vc, map, std::make_pair(ibeg+offset-nx, iend+offset-nx),
                            std::make_pair(jbeg   , jend),
                            std::make_pair(kbeg   , kend));

  // Load face-centered field in the buffer
  if(haveVs) {
    BufferLeft.Pack(Vs, BX1s,std::make_pair(ibeg+nx+1, iend+nx+1),
                             std::make_pair(jbeg   , jend),
                             std::make_pair(kbeg   , kend));

    BufferRight.Pack(Vs, BX1s, std::make_pair(ibeg+offset-nx, iend+offset-nx),
                               std::make_pair(jbeg   , jend),
                               std::make_pair(kbeg   , kend));

    #if DIMENSIONS >= 2

    BufferLeft.Pack(Vs, BX2s,std::make_pair(ibeg+nx, iend+nx),
                             std::make_pair(jbeg   , jend+1),
                             std::make_pair(kbeg   , kend));

    BufferRight.Pack(Vs, BX2s, std::make_pair(ibeg+offset-nx, iend+offset-nx),
                               std::make_pair(jbeg   , jend+1),
                               std::make_pair(kbeg   , kend));

    #endif

    #if DIMENSIONS == 3

    BufferLeft.Pack(Vs, BX3s,std::make_pair(ibeg+nx, iend+nx),
                             std::make_pair(jbeg   , jend),
                             std::make_pair(kbeg   , kend+1));

    BufferRight.Pack(Vs, BX3s, std::make_pair(ibeg+offset-nx, iend+offset-nx),
                               std::make_pair(jbeg   , jend),
                               std::make_pair(kbeg   , kend+1));

    #endif
  }

  // Wait for completion before sending out everything
  Kokkos::fence();
  myTimer -= MPI_Wtime();
  tStart = MPI_Wtime();
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
  myTimer += MPI_Wtime();
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;
  // Unpack
  BufferLeft=BufferRecvX1[faceLeft];
  BufferRight=BufferRecvX1[faceRight];

  BufferLeft.ResetPointer();
  BufferRight.ResetPointer();

  BufferLeft.Unpack(Vc, map,std::make_pair(ibeg, iend),
                            std::make_pair(jbeg   , jend),
                            std::make_pair(kbeg   , kend));

  BufferRight.Unpack(Vc, map,std::make_pair(ibeg+offset, iend+offset),
                             std::make_pair(jbeg   , jend),
                             std::make_pair(kbeg   , kend));
  // We fill the ghost zones

  if(haveVs) {
    BufferLeft.Unpack(Vs, BX1s, std::make_pair(ibeg, iend),
                                std::make_pair(jbeg   , jend),
                                std::make_pair(kbeg   , kend));

    BufferRight.Unpack(Vs, BX1s, std::make_pair(ibeg+offset+1, iend+offset+1),
                                std::make_pair(jbeg   , jend),
                                std::make_pair(kbeg   , kend));

    #if DIMENSIONS >= 2
    BufferLeft.Unpack(Vs, BX2s, std::make_pair(ibeg, iend),
                                std::make_pair(jbeg   , jend+1),
                                std::make_pair(kbeg   , kend));

    BufferRight.Unpack(Vs, BX2s, std::make_pair(ibeg+offset, iend+offset),
                                std::make_pair(jbeg   , jend+1),
                                std::make_pair(kbeg   , kend));
    #endif

    #if DIMENSIONS == 3
    BufferLeft.Unpack(Vs, BX3s, std::make_pair(ibeg, iend),
                                std::make_pair(jbeg   , jend),
                                std::make_pair(kbeg   , kend+1));

    BufferRight.Unpack(Vs, BX3s, std::make_pair(ibeg+offset, iend+offset),
                                std::make_pair(jbeg   , jend),
                                std::make_pair(kbeg   , kend+1));
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
  int ibeg,iend,jbeg,jend,kbeg,kend,offset,ny;
  Buffer BufferLeft=BufferSendX2[faceLeft];
  Buffer BufferRight=BufferSendX2[faceRight];
  IdefixArray1D<int> map = this->mapVars;

// If MPI Persistent, start receiving even before the buffers are filled
  myTimer -= MPI_Wtime();
  double tStart = MPI_Wtime();
#ifdef MPI_PERSISTENT
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];

  MPI_SAFE_CALL(MPI_Startall(2, recvRequestX2));
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;
#endif
  myTimer += MPI_Wtime();

  // Coordinates of the ghost region which needs to be transfered
  ibeg   = 0;
  iend   = ntot[IDIR];

  jbeg   = 0;
  jend   = nghost[JDIR];
  offset = end[JDIR];     // Distance between beginning of left and right ghosts
  ny     = nghost[JDIR];

  kbeg   = beg[KDIR];
  kend   = end[KDIR];

  BufferLeft.ResetPointer();
  BufferRight.ResetPointer();

  BufferLeft.Pack(Vc, map, std::make_pair(ibeg    , iend),
                           std::make_pair(jbeg+ny , jend+ny),
                           std::make_pair(kbeg    , kend));

  BufferRight.Pack(Vc, map, std::make_pair(ibeg           , iend),
                            std::make_pair(jbeg+offset-ny , jend+offset-ny),
                            std::make_pair(kbeg           , kend));

  // Load face-centered field in the buffer
  if(haveVs) {
    BufferLeft.Pack(Vs, BX1s,std::make_pair(ibeg , iend+1),
                          std::make_pair(jbeg+ny , jend+ny),
                          std::make_pair(kbeg    , kend));

    BufferRight.Pack(Vs, BX1s, std::make_pair(ibeg        , iend+1),
                            std::make_pair(jbeg+offset-ny , jend+offset-ny),
                            std::make_pair(kbeg           , kend));
    #if DIMENSIONS >= 2
    BufferLeft.Pack(Vs, BX2s,std::make_pair(ibeg , iend),
                             std::make_pair(jbeg+ny+1 , jend+ny+1),
                             std::make_pair(kbeg    , kend));

    BufferRight.Pack(Vs, BX2s, std::make_pair(ibeg        , iend),
                            std::make_pair(jbeg+offset-ny , jend+offset-ny),
                            std::make_pair(kbeg           , kend));
    #endif
    #if DIMENSIONS == 3

    BufferLeft.Pack(Vs, BX3s,std::make_pair(ibeg , iend),
                          std::make_pair(jbeg+ny , jend+ny),
                          std::make_pair(kbeg    , kend+1));

    BufferRight.Pack(Vs, BX3s, std::make_pair(ibeg        , iend),
                            std::make_pair(jbeg+offset-ny , jend+offset-ny),
                            std::make_pair(kbeg           , kend+1));

    #endif
  }

  // Send to the right
  Kokkos::fence();

  myTimer -= MPI_Wtime();
  tStart = MPI_Wtime();
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
  myTimer += MPI_Wtime();
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;
  // Unpack
  BufferLeft=BufferRecvX2[faceLeft];
  BufferRight=BufferRecvX2[faceRight];

  BufferLeft.ResetPointer();
  BufferRight.ResetPointer();

  // We fill the ghost zones
  BufferLeft.Unpack(Vc, map,std::make_pair(ibeg, iend),
                            std::make_pair(jbeg   , jend),
                            std::make_pair(kbeg   , kend));

  BufferRight.Unpack(Vc, map,std::make_pair(ibeg        , iend),
                             std::make_pair(jbeg+offset , jend+offset),
                             std::make_pair(kbeg        , kend));
  // We fill the ghost zones

  if(haveVs) {
    BufferLeft.Unpack(Vs, BX1s, std::make_pair(ibeg, iend+1),
                                std::make_pair(jbeg   , jend),
                                std::make_pair(kbeg   , kend));

    BufferRight.Unpack(Vs, BX1s, std::make_pair(ibeg, iend+1),
                                std::make_pair(jbeg+offset   , jend+offset),
                                std::make_pair(kbeg   , kend));
    #if DIMENSIONS >= 2
    BufferLeft.Unpack(Vs, BX2s, std::make_pair(ibeg, iend),
                                std::make_pair(jbeg   , jend),
                                std::make_pair(kbeg   , kend));

    BufferRight.Unpack(Vs, BX2s, std::make_pair(ibeg, iend),
                                std::make_pair(jbeg+offset+1, jend+offset+1),
                                std::make_pair(kbeg   , kend));
    #endif
    #if DIMENSIONS == 3
    BufferLeft.Unpack(Vs, BX3s, std::make_pair(ibeg, iend),
                                std::make_pair(jbeg, jend),
                                std::make_pair(kbeg, kend+1));

    BufferRight.Unpack(Vs, BX3s,std::make_pair(ibeg      , iend),
                                std::make_pair(jbeg+offset, jend+offset),
                                std::make_pair(kbeg       , kend+1));
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
  int ibeg,iend,jbeg,jend,kbeg,kend,offset,nz;
  Buffer BufferLeft=BufferSendX3[faceLeft];
  Buffer BufferRight=BufferSendX3[faceRight];
  IdefixArray1D<int> map = this->mapVars;

  // If MPI Persistent, start receiving even before the buffers are filled
  myTimer -= MPI_Wtime();

  double tStart = MPI_Wtime();
#ifdef MPI_PERSISTENT
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];

  MPI_SAFE_CALL(MPI_Startall(2, recvRequestX3));
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;
#endif
  myTimer += MPI_Wtime();
  // Coordinates of the ghost region which needs to be transfered
  ibeg   = 0;
  iend   = ntot[IDIR];

  jbeg   = 0;
  jend   = ntot[JDIR];

  kbeg   = 0;
  kend   = nghost[KDIR];
  offset = end[KDIR];     // Distance between beginning of left and right ghosts
  nz     = nghost[KDIR];

  BufferLeft.ResetPointer();
  BufferRight.ResetPointer();

  BufferLeft.Pack(Vc, map, std::make_pair(ibeg   , iend),
                           std::make_pair(jbeg   , jend),
                           std::make_pair(kbeg+nz, kend+nz));

  BufferRight.Pack(Vc, map, std::make_pair(ibeg            , iend),
                            std::make_pair(jbeg            , jend),
                            std::make_pair(kbeg + offset-nz, kend+ offset-nz));

  // Load face-centered field in the buffer
  if(haveVs) {
    BufferLeft.Pack(Vs, BX1s,std::make_pair(ibeg , iend+1),
                             std::make_pair(jbeg    , jend),
                             std::make_pair(kbeg+nz , kend+nz));

    BufferRight.Pack(Vs, BX1s, std::make_pair(ibeg            , iend+1),
                               std::make_pair(jbeg            , jend),
                               std::make_pair(kbeg + offset-nz, kend+ offset-nz));

    #if DIMENSIONS >= 2

    BufferLeft.Pack(Vs, BX2s,std::make_pair(ibeg    , iend),
                             std::make_pair(jbeg    , jend+1),
                             std::make_pair(kbeg+nz , kend+nz));

    BufferRight.Pack(Vs, BX2s, std::make_pair(ibeg            , iend),
                               std::make_pair(jbeg            , jend+1),
                               std::make_pair(kbeg + offset-nz, kend+ offset-nz));

    #endif

    #if DIMENSIONS == 3
    BufferLeft.Pack(Vs, BX3s,std::make_pair(ibeg    , iend),
                             std::make_pair(jbeg    , jend),
                             std::make_pair(kbeg+nz+1 , kend+nz+1));

    BufferRight.Pack(Vs, BX3s, std::make_pair(ibeg            , iend),
                               std::make_pair(jbeg            , jend),
                               std::make_pair(kbeg + offset-nz, kend+ offset-nz));
    #endif
  }

  // Send to the right
  Kokkos::fence();

  myTimer -= MPI_Wtime();
  tStart = MPI_Wtime();
#ifdef MPI_PERSISTENT
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

  MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX3[faceLeft].data(), bufferSizeX3, realMPI, procSend, 301,
                BufferRecvX3[faceRight].data(), bufferSizeX3, realMPI, procRecv, 301,
                mygrid->CartComm, &status));
  #endif
#endif
  myTimer += MPI_Wtime();
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;
  // Unpack
  BufferLeft=BufferRecvX3[faceLeft];
  BufferRight=BufferRecvX3[faceRight];

  BufferLeft.ResetPointer();
  BufferRight.ResetPointer();


  // We fill the ghost zones
  BufferLeft.Unpack(Vc, map,std::make_pair(ibeg, iend),
                            std::make_pair(jbeg   , jend),
                            std::make_pair(kbeg   , kend));

  BufferRight.Unpack(Vc, map,std::make_pair(ibeg        , iend),
                             std::make_pair(jbeg        , jend),
                             std::make_pair(kbeg+offset , kend+offset));
  // We fill the ghost zones

  if(haveVs) {
    BufferLeft.Unpack(Vs, BX1s, std::make_pair(ibeg, iend+1),
                                std::make_pair(jbeg   , jend),
                                std::make_pair(kbeg   , kend));

    BufferRight.Unpack(Vs, BX1s, std::make_pair(ibeg, iend+1),
                                std::make_pair(jbeg          , jend),
                                std::make_pair(kbeg+offset   , kend+offset));

    #if DIMENSIONS >=2
    BufferLeft.Unpack(Vs, BX2s, std::make_pair(ibeg, iend),
                                std::make_pair(jbeg, jend+1),
                                std::make_pair(kbeg, kend));

    BufferRight.Unpack(Vs, BX2s,std::make_pair(ibeg       , iend),
                                std::make_pair(jbeg       , jend+1),
                                std::make_pair(kbeg+offset, kend+offset));
    #endif

    #if DIMENSIONS == 3
    BufferLeft.Unpack(Vs, BX3s, std::make_pair(ibeg, iend),
                                std::make_pair(jbeg, jend),
                                std::make_pair(kbeg, kend));

    BufferRight.Unpack(Vs, BX3s,std::make_pair(ibeg       , iend),
                                std::make_pair(jbeg       , jend),
                                std::make_pair(kbeg+offset+1, kend+offset+1));
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
  bytesSentOrReceived += 4*bufferSizeX3*sizeof(real);

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

  if(idfx::prank == 0) {
    srcHost(0) = 0;
    Kokkos::deep_copy(src, srcHost);
  }

  if(idfx::psize > 1) {
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

    // Capture segfaults
    struct sigaction newHandler;
    struct sigaction oldHandler;
    memset(&newHandler, 0, sizeof(newHandler));
    newHandler.sa_flags = SA_SIGINFO;
    newHandler.sa_sigaction = Mpi::SigErrorHandler;
    sigaction(SIGSEGV, &newHandler, &oldHandler);
    try {
      // We next circulate the info round-robin accross all the nodes to check that
      // MPI can exchange buffers in idefix arrays

      MPI_Status status;
      int ierrSend, ierrRecv;
      if(idfx::prank == 0) {
        ierrSend = MPI_Send(src.data(), 1, MPI_INT64_T, idfx::prank+1, 1, MPI_COMM_WORLD);
        ierrRecv = MPI_Recv(src.data(), 1, MPI_INT64_T, idfx::psize-1, 1, MPI_COMM_WORLD, &status);
      } else {
        ierrRecv = MPI_Recv(src.data(), 1, MPI_INT64_T, idfx::prank-1, 1, MPI_COMM_WORLD, &status);
        // Add our own rank to the data
        Kokkos::deep_copy(srcHost, src);
        srcHost(0) += idfx::prank;
        Kokkos::deep_copy(src, srcHost);
        ierrSend = MPI_Send(src.data(), 1, MPI_INT64_T, (idfx::prank+1)%idfx::psize, 1,
                                                        MPI_COMM_WORLD);
      }

      if(ierrSend != 0) {
        char MPImsg[MPI_MAX_ERROR_STRING];
        int MPImsgLen;
        MPI_Error_string(ierrSend, MPImsg, &MPImsgLen);
        throw std::runtime_error(std::string(MPImsg, MPImsgLen));
      }
      if(ierrRecv != 0) {
        char MPImsg[MPI_MAX_ERROR_STRING];
        int MPImsgLen;
        MPI_Error_string(ierrSend, MPImsg, &MPImsgLen);
        throw std::runtime_error(std::string(MPImsg, MPImsgLen));
      }
    } catch(std::exception &e) {
      std::stringstream errmsg;
      errmsg << "Your MPI library is unable to perform Send/Recv on Idefix arrays.";
      errmsg << std::endl;
      #ifdef KOKKOS_ENABLE_CUDA
        errmsg << "Check that your MPI library is CUDA aware." << std::endl;
      #elif defined(KOKKOS_ENABLE_HIP)
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
  }

  // Check that we have the proper end result
  Kokkos::deep_copy(srcHost, src);
  int64_t size = static_cast<int64_t>(idfx::psize);
  int64_t rank = static_cast<int64_t>(idfx::prank);
  int64_t result = rank == 0 ? size*(size-1)/2 : rank*(rank+1)/2;

  if(srcHost(0) != result) {
    idfx::cout << "got " << srcHost(0) << " expected " << result << std::endl;
    std::stringstream errmsg;
    errmsg << "Your MPI library managed to perform MPI exchanges on Idefix Arrays, but the result ";
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
  #elif defined(KOKKOS_ENABLE_HIP)
    errmsg << "Check that your MPI library is RocM aware." << std::endl;
  #else
    errmsg << "Check your MPI library configuration." << std::endl;
  #endif
  IDEFIX_ERROR(errmsg);
}

// This routine check that all of the processes are synced.
// Returns true if this is the case, false otherwise

bool Mpi::CheckSync(real timeout) {
  // If no parallelism, then we're in sync!
  if(idfx::psize == 1) return(true);

  int send = idfx::prank;
  int recv = 0;
  MPI_Request request;

  MPI_Iallreduce(&send, &recv, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, &request);

  double start = MPI_Wtime();
  int flag = 0;
  MPI_Status status;

  while((MPI_Wtime()-start < timeout) && !flag) {
    MPI_Test(&request, &flag, &status);
    // sleep for 10 ms
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }
  if(!flag) {
    // We did not managed to do an allreduce, so this is a failure.
    return(false);
  }
  if(recv != idfx::psize*(idfx::psize-1)/2) {
    IDEFIX_ERROR("wrong result for synchronisation");
  }

  return(true);
}
