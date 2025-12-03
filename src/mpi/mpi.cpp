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

  // increase the number of instances
  nInstances++;
  thisInstance=nInstances;

  for(int dir=0; dir<3; dir++) {
    std::array<bool,2> overWriteBXn = {true, true};
    if(grid->lbound[dir] == BoundaryType::shearingbox) {
      overWriteBXn[faceLeft] = false;
    }
    if(grid->rbound[dir] == BoundaryType::shearingbox) {
      overWriteBXn[faceRight] = false;
    }

    exchanger[dir].Init(grid, dir, inputMap,
                          {nghost[0], nghost[1], nghost[2]},
                          {nint[0], nint[1], nint[2]},
                          inputHaveVs, overWriteBXn);
  }

  isInitialized = true;
  idfx::popRegion();
}

// Destructor (clean up persistent communication channels)
Mpi::~Mpi() {
  idfx::pushRegion("Mpi::~Mpi");
  if(isInitialized) {
    if(thisInstance==1) {
      int bytesSentOrReceived = 0;
      double myTimer = 0;
      for(int dir=0; dir<3; dir++) {
        bytesSentOrReceived += exchanger[dir].bytesSentOrReceived;
        myTimer += exchanger[dir].myTimer;
      }
      idfx::cout << "Mpi(" << thisInstance << "): measured throughput is "
                << bytesSentOrReceived/myTimer/1024.0/1024.0 << " MB/s" << std::endl;
      idfx::cout << "Mpi(" << thisInstance << "): message sizes were " << std::endl;
      idfx::cout << "        X1: "
                 << exchanger[IDIR].bufferSize[0]*sizeof(real)/1024.0/1024.0
                 << " MB" << std::endl;
      idfx::cout << "        X2: "
                 << exchanger[JDIR].bufferSize[0]*sizeof(real)/1024.0/1024.0
                 << " MB" << std::endl;
      idfx::cout << "        X3: "
                 << exchanger[KDIR].bufferSize[0]*sizeof(real)/1024.0/1024.0
                 << " MB" << std::endl;
    }
    isInitialized = false;
  }
  idfx::popRegion();
}

void Mpi::ExchangeX1(IdefixArray4D<real> Vc, IdefixArray4D<real> Vs) {
  idfx::pushRegion("Mpi::ExchangeX1");

  exchanger[IDIR].Exchange(Vc, Vs);
  idfx::popRegion();
}

void Mpi::ExchangeX2(IdefixArray4D<real> Vc, IdefixArray4D<real> Vs) {
  idfx::pushRegion("Mpi::ExchangeX2");
  exchanger[JDIR].Exchange(Vc, Vs);
  idfx::popRegion();
}

void Mpi::ExchangeX3(IdefixArray4D<real> Vc, IdefixArray4D<real> Vs) {
  idfx::pushRegion("Mpi::ExchangeX3");
  exchanger[KDIR].Exchange(Vc, Vs);
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
