// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <vector>
#include "exchanger.hpp"
#include "idefix.hpp"
#include "buffer.hpp"
#include "grid.hpp"
#include "arrays.hpp"

//#define MPI_NON_BLOCKING
#define MPI_PERSISTENT

int Exchanger::nInstances = 0;

void Exchanger::Init(
              Grid *grid,
              int direction,
              std::vector<int> inputMap,
              std::array<int, 3> nghost,
              std::array<int, 3> nint,
              bool inputHaveVs,
              std::array<bool,2> overwriteBXn) {
  this->grid = grid;
  this->direction = direction;
  // Allocate mapVars on target and copy it from the input argument list
  this->mapVars = idfx::ConvertVectorToIdefixArray(inputMap);
  this->mapNVars = inputMap.size();
  this->haveVs = inputHaveVs;

  // increase the number of instances
  this->thisInstance = nInstances;

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
  // Buffer size direction are for sends (i.e. buffer left is for a send from the left side)

  // Make left buffer
  // Init zone to the full domain
  for(int dir = 0 ; dir < 3 ; dir++) {
    if(dir < direction) {
      boxRecv[faceLeft][dir][0] = 0;
      boxRecv[faceLeft][dir][1] = ntot[dir];
    } else if(dir > direction) {
      boxRecv[faceLeft][dir][0] = beg[dir];
      boxRecv[faceLeft][dir][1] = end[dir];
    } else {
      // dir == direction
      boxRecv[faceLeft][dir][0] = 0;
      boxRecv[faceLeft][dir][1] = nghost[dir];
    }
  }
  // Copy all the boxes from boxRecvLeft
  boxRecv[faceRight] = boxRecv[faceLeft];
  boxSend = boxRecv;

  // Adjust the indices for send and receive in the direction of interest
  boxRecv[faceRight][direction][0] = end[direction];
  boxRecv[faceRight][direction][1] = ntot[direction];

  boxSend[faceLeft][direction][0] = beg[direction];
  boxSend[faceLeft][direction][1] = beg[direction]+nghost[direction];

  boxSend[faceRight][direction][0] = end[direction] - nghost[direction];
  boxSend[faceRight][direction][1] = end[direction];

  // Face-centered boxes

  // Add one element in the field direction
  for(int component = 0 ; component < 3 ; component++) {
    // Init as centered boxes
    boxSendVs[component] = boxSend;
    boxRecvVs[component] = boxRecv;
    const int normalDir = component;
    if(component != direction) {
      for(int face = 0 ; face <2 ; face++) {
        boxSendVs[component][face][normalDir][1] += 1;
        boxRecvVs[component][face][normalDir][1] += 1;
      }
    } else {
      // component == direction
      if(!overwriteBXn[faceLeft]) boxSendVs[component][faceLeft][normalDir][0] += 1;
      boxSendVs[component][faceLeft][normalDir][1] += 1;

      if(!overwriteBXn[faceRight]) boxRecvVs[component][faceRight][normalDir][0] += 1;
      boxRecvVs[component][faceRight][normalDir][1] += 1;
    }
  }

  // Compute buffer sizes
  for(int face=0 ; face < 2 ; face++) {
    bufferSizeSend[face] = mapNVars * Buffer::ComputeBoxSize(boxSend[face]);
    bufferSizeRecv[face] = mapNVars * Buffer::ComputeBoxSize(boxRecv[face]);
    if(haveVs) {
      for(int component = 0 ; component <DIMENSIONS ; component++) {
        bufferSizeSend[face] += Buffer::ComputeBoxSize( boxSendVs[component][face] );
        bufferSizeRecv[face] += Buffer::ComputeBoxSize( boxRecvVs[component][face] );
      }
    }
  }

  // allocate buffers
  BufferSend[faceLeft] = Buffer(bufferSizeSend[faceLeft]);
  BufferRecv[faceRight] = Buffer(bufferSizeRecv[faceRight]);

  BufferSend[faceRight] = Buffer(bufferSizeSend[faceRight]);
  BufferRecv[faceLeft] = Buffer(bufferSizeRecv[faceLeft]);

  // Compute directions
  MPI_Cart_shift(grid->CartComm,direction,1,&procRecv[faceLeft],&procSend[faceRight]);
  MPI_Cart_shift(grid->CartComm,direction,-1,&procRecv[faceRight],&procSend[faceLeft]);

  #ifdef MPI_PERSISTENT

  // X1-dir exchanges
  // We receive from procRecv, and we send to procSend

  MPI_Send_init(BufferSend[faceRight].data(), bufferSizeSend[faceRight], realMPI,
            procSend[faceRight], thisInstance*2,
            grid->CartComm, &sendRequest[faceRight]);

  MPI_Recv_init(BufferRecv[faceLeft].data(), bufferSizeRecv[faceLeft], realMPI,
            procRecv[faceLeft],thisInstance*2,
            grid->CartComm, &recvRequest[faceLeft]);

  // Send to the left
  // We receive from procRecv, and we send to procSend

  MPI_Send_init(BufferSend[faceLeft].data(), bufferSizeSend[faceLeft], realMPI,
            procSend[faceLeft],thisInstance*2+1,
            grid->CartComm, &sendRequest[faceLeft]);

  MPI_Recv_init(BufferRecv[faceRight].data(), bufferSizeRecv[faceRight], realMPI,
            procRecv[faceRight], thisInstance*2+1,
            grid->CartComm, &recvRequest[faceRight]);

  #endif // MPI_PERSISTENT

  // say this instance is initialized.
  isInitialized = true;
  nInstances++;

  idfx::popRegion();
}

Exchanger::~Exchanger() {
  idfx::pushRegion("Exchanger::~Exchanger");
  if(isInitialized) {
    // Properly clean up the mess
    #ifdef MPI_PERSISTENT
      for(int i=0 ; i< 2; i++) {
        MPI_Request_free( &sendRequest[i]);
        MPI_Request_free( &recvRequest[i]);
      }
    #endif
    isInitialized = false;
  }
  idfx::popRegion();
}

void Exchanger::Exchange(IdefixArray4D<real> Vc, IdefixArray4D<real> Vs) {
  idfx::pushRegion("Mpi::ExchangeX1");
  // Load  the buffers with data
  Buffer BufferLeft = BufferSend[faceLeft];
  Buffer BufferRight = BufferSend[faceRight];
  IdefixArray1D<int> map = this->mapVars;

  bool recvRight = (procRecv[faceRight] != MPI_PROC_NULL);
  bool recvLeft  = (procRecv[faceLeft] != MPI_PROC_NULL);

  // If MPI Persistent, start receiving even before the buffers are filled
  myTimer -= MPI_Wtime();
  double tStart = MPI_Wtime();
#ifdef MPI_PERSISTENT
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];

  MPI_Startall(2, recvRequest);
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;
#endif
  myTimer += MPI_Wtime();

  BufferLeft.ResetPointer();
  BufferRight.ResetPointer();

  BufferLeft.Pack(Vc, map, boxSend[faceLeft]);
  BufferRight.Pack(Vc, map, boxSend[faceRight]);
  // Load face-centered field in the buffer
  if(haveVs) {
    for(int component = 0 ; component < DIMENSIONS ; component++) {
      BufferLeft.Pack(Vs, component, boxSendVs[component][faceLeft]);
      BufferRight.Pack(Vs, component, boxSendVs[component][faceRight]);
    }
  }

  // Wait for completion before sending out everything
  Kokkos::fence();
  myTimer -= MPI_Wtime();
  tStart = MPI_Wtime();
#ifdef MPI_PERSISTENT
  MPI_Startall(2, sendRequest);
  // Wait for buffers to be received
  MPI_Waitall(2,recvRequest,recvStatus);

#else

  #ifdef MPI_NON_BLOCKING
  MPI_Status sendStatus[2];
  MPI_Status recvStatus[2];
  MPI_Request sendRequest[2];
  MPI_Request recvRequest[2];

  // We receive from procRecv, and we send to procSend

  MPI_Isend(BufferSend[faceRight].data(), bufferSizeSend[faceRight], realMPI,
            procSend[faceRight], 100, mygrid->CartComm, &sendRequest[0]);

  MPI_Irecv(BufferRecv[faceLeft].data(), bufferSizeRecv[faceLeft], realMPI,
            procRecv[faceLeft], 100, mygrid->CartComm, &recvRequest[0]);
  // Send to the left
  // We receive from procRecv, and we send to procSend

  MPI_Isend(BufferSend[faceLeft].data(), bufferSizeSend[faceLeft], realMPI,
            procSend[faceLeft], 101, mygrid->CartComm, &sendRequest[1]);

  MPI_Irecv(BufferRecv[faceRight].data(), bufferSizeRecv[faceRight], realMPI,
            procRecv[faceRight], 101, mygrid->CartComm, &recvRequest[1]);

  // Wait for recv to complete (we don't care about the sends)
  MPI_Waitall(2, recvRequest, recvStatus);

  #else
  MPI_Status status;
  // Send to the right
  // We receive from procRecv, and we send to procSend

  MPI_Sendrecv(BufferSend[faceRight].data(), bufferSizeSend[faceRight], realMPI,
                procSend[faceRight], 100,
                BufferRecv[faceLeft].data(), bufferSizeRecv[faceLeft], realMPI,
                procRecv[faceLeft], 100,
                grid->CartComm, &status);

  // Send to the left
  // We receive from procRecv, and we send to procSend

  MPI_Sendrecv(BufferSend[faceLeft].data(), bufferSizeSend[faceLeft], realMPI,
                procSend[faceLeft], 101,
                BufferRecv[faceRight].data(), bufferSizeRecv[faceRight], realMPI,
                procRecv[faceRight], 101,
                grid->CartComm, &status);
  #endif
#endif
myTimer += MPI_Wtime();
idfx::mpiCallsTimer += MPI_Wtime() - tStart;
// Unpack
BufferLeft=BufferRecv[faceLeft];
BufferRight=BufferRecv[faceRight];

BufferLeft.ResetPointer();
BufferRight.ResetPointer();

if(recvLeft) {
    BufferLeft.Unpack(Vc, map, boxRecv[faceLeft]);
  if(haveVs) {
    for(int component = 0 ; component < DIMENSIONS ; component++) {
      BufferLeft.Unpack(Vs, component, boxRecvVs[component][faceLeft]);
    }
  }
}
if(recvRight) {
  BufferRight.Unpack(Vc, map, boxRecv[faceRight]);
  if(haveVs) {
    for(int component = 0 ; component < DIMENSIONS ; component++) {
      BufferRight.Unpack(Vs, component, boxRecvVs[component][faceRight]);
    }
  }
}
myTimer -= MPI_Wtime();
#ifdef MPI_NON_BLOCKING
  // Wait for the sends if they have not yet completed
  MPI_Waitall(2, sendRequest, sendStatus);
#endif

#ifdef MPI_PERSISTENT
  MPI_Waitall(2, sendRequest, sendStatus);
#endif
  myTimer += MPI_Wtime();
  bytesSentOrReceived += (bufferSizeRecv[faceLeft]
                          +bufferSizeSend[faceLeft]
                          +bufferSizeRecv[faceRight]
                          +bufferSizeSend[faceRight])*sizeof(real);

  idfx::popRegion();
}
