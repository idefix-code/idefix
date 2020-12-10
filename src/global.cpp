// ***********************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************************

#include <string>
#include <sstream>
#include "idefix.hpp"
#include "global.hpp"

namespace idfx {

int prank;
int psize;

double mpiTimer;

IdefixOstream cout;

#ifdef DEBUG
static int regionIndent = 0;
#endif

#ifdef WITH_MPI
extern MPI_Comm CartComm;
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
  mpiTimer = 0.0;
  return(0);
}   // Initialisation routine for idefix

void pushRegion(const std::string& kName) {
  Kokkos::Profiling::pushRegion(kName);

#ifdef DEBUG
  regionIndent=regionIndent+4;
  for(int i=0; i < regionIndent ; i++) {
    cout << "-";
  }
  cout << "> " << kName << std::endl;
#endif
}

void popRegion() {
  Kokkos::Profiling::popRegion();
#ifdef DEBUG
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

} // namespace idfx
