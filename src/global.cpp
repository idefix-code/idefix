// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <iostream>
#include <string>
#include <sstream>
#include "idefix.hpp"
#include "global.hpp"
#include "profiler.hpp"

#ifdef WITH_MPI
#include "mpi.hpp"
#endif

namespace idfx {

int prank;
int psize;

double mpiCallsTimer = 0.0;

bool warningsAreErrors{false};

IdefixOutStream cout;
IdefixErrStream cerr;
Profiler prof;
LoopPattern defaultLoopPattern;

#ifdef DEBUG
static int regionIndent = 0;
#endif

int initialize() {
#ifdef WITH_MPI
  MPI_Comm_size(MPI_COMM_WORLD,&psize);
  MPI_Comm_rank(MPI_COMM_WORLD,&prank);
#else
  psize=1;
  prank=0;
#endif
  cout.init(prank);
  prof.Init();

  // Init loop Pattern
  #if defined(KOKKOS_ENABLE_CUDA)
    defaultLoopPattern = LoopPattern::RANGE;  // On cuda, works best (generally)
  #elif defined(KOKKOS_ENABLE_HIP)
    defaultLoopPattern = LoopPattern::RANGE;  // On HIP, works best (generally)
  #else
    defaultLoopPattern = LoopPattern::TPX;    // On cpus, works best (generally)
  #endif

  #ifdef WITH_MPI
    Mpi::CheckConfig();
  #endif

  return(0);
}   // Initialisation routine for idefix

void pushRegion(const std::string& kName) {
  Kokkos::Profiling::pushRegion(kName);
  if(prof.perfEnabled) {
    prof.currentRegion = prof.currentRegion->GetChild(kName);
    prof.currentRegion->Start();
  }
#ifdef DEBUG
  regionIndent=regionIndent+4;
  for(int i=0; i < regionIndent ; i++) {
    cout << "-";
  }
  cout << "> " << kName << "..." << std::endl;
#endif
}

void popRegion() {
  Kokkos::Profiling::popRegion();
  if(prof.perfEnabled) {
    Kokkos::fence();
    prof.currentRegion->Stop();
    prof.currentRegion = prof.currentRegion->parent;
  }
#ifdef DEBUG
  for(int i=0; i < regionIndent ; i++) {
    cout << "-";
  }
  cout << "> ...returned" << std::endl;
  regionIndent = regionIndent-4;
#endif
}

// Init the iostream with defined rank
void IdefixOutStream::init(int rank) {
  if(rank==0)
    this->toscreen=true;
  else
    this->toscreen=false;
}


// disable the log file
void IdefixOutStream::enableLogFile() {
  std::stringstream sslogFileName;
  sslogFileName << "idefix." << idfx::prank << ".log";

  std::string logFileName(sslogFileName.str());
  this->my_fstream.open(logFileName.c_str());

  this->logFileEnabled = true;
}


/*********************************************/
/**
Customized random number generator
Allow one to have consistant random numbers
generators on different architectures.
**/
/*********************************************/
real randm(void) {
    const int a = 16807;
    const int m = 2147483647;
    static int in0 = 13763 + 2417*prank;
    int q;

    /* find random number  */
    q = static_cast<int>(fmod(static_cast<double>(a) * in0, m));
    in0 = q;

    return static_cast<real>(static_cast<double>(q) / static_cast<double>(m));
}

void safeExit(int retCode) {
  if(retCode != 0) {
    #ifdef WITH_MPI
    MPI_Abort(MPI_COMM_WORLD,retCode);
    #endif
  } else {
    Kokkos::finalize();
    #ifdef WITH_MPI
    MPI_Finalize();
    #endif
  }
  exit(retCode);
}

} // namespace idfx
