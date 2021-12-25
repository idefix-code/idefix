// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>
#include <sstream>
#include "idefix.hpp"
#include "global.hpp"
#include "profiler.hpp"

namespace idfx {

int prank;
int psize;

double mpiCallsTimer = 0.0;

bool warningsAreErrors{false};

IdefixOstream cout;
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

  return(0);
}   // Initialisation routine for idefix

void pushRegion(const std::string& kName) {
  Kokkos::Profiling::pushRegion(kName);

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
#ifdef DEBUG
  for(int i=0; i < regionIndent ; i++) {
    cout << "-";
  }
  cout << "> ...returned" << std::endl;
  regionIndent = regionIndent-4;
#endif
}

// Init the iostream with defined rank
void IdefixOstream::init(int rank) {
  std::stringstream sslogFileName;
  sslogFileName << "idefix." << rank << ".log";

  std::string logFileName(sslogFileName.str());
  this->my_fstream.open(logFileName.c_str());

  if(rank==0)
    this->toscreen=true;
  else
    this->toscreen=false;
}

// disable the log file
void IdefixOstream::disableLogFile() {
  my_fstream << "Log files have been disabled (e.g. using -nowrite)." << std::endl;
  my_fstream.close();
  this->logFileEnabled = false;
}

} // namespace idfx
