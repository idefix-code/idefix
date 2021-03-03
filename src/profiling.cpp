// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <atomic>
#include <mutex>
#include "idefix.hpp"
#include "profiling.hpp"

void idfx::Profiler::Init() {
  idfx::pushRegion("Profiler::Init");
  idfx::cout << "init" << std::endl;
  this->numSpaces=0;
  for(int i=0; i < 16 ; i++) {
    this->spaceSize[i] = 0;
    this->spaceMax[i] = 0;
  }
  idfx::popRegion();
}

void idfx::Profiler::Show() {
  idfx::pushRegion("Profiler::Show");
  idfx::cout << "coucou" << std::endl;
  this->numSpaces=0;
  for(int i=0; i < this->numSpaces ; i++) {
    idfx::cout << "Profiler: maximum memory usage for space " << this->spaceName[i]
               << " : " << this->spaceMax[i]/(1024*1024) << "MB." << std::endl;
  }
  idfx::popRegion();
}


// Kokkos Profiler hooks

extern "C" void kokkosp_init_library(const int loadSeq,
  const uint64_t interfaceVer,
  const uint32_t devInfoCount,
  void* deviceInfo) {

  // Init profiling class
  idfx::prof.Init();
}

extern "C" void kokkosp_allocate_data(const idfx::SpaceHandle space, const char* label, const void* const ptr, const uint64_t size) {
  std::lock_guard<std::mutex> lock(idfx::prof.m);

  int space_i = idfx::prof.numSpaces;
  for(int s = 0; s<idfx::prof.numSpaces; s++)
    if(strcmp(idfx::prof.spaceName[s],space.name)==0)
      space_i = s;

  if(space_i == idfx::prof.numSpaces) {
    strncpy(idfx::prof.spaceName[space_i],space.name,64);
    idfx::prof.numSpaces++;
  }
  idfx::prof.spaceSize[space_i] += size;
  if(idfx::prof.spaceSize[space_i] > idfx::prof.spaceMax[space_i]) {
    idfx::prof.spaceMax[space_i] = idfx::prof.spaceSize[space_i];
  }
}


extern "C" void kokkosp_deallocate_data(const idfx::SpaceHandle space, const char* label, const void* const ptr, const uint64_t size) {
  std::lock_guard<std::mutex> lock(idfx::prof.m);

  int space_i = idfx::prof.numSpaces;
  for(int s = 0; s<idfx::prof.numSpaces; s++)
    if(strcmp(idfx::prof.spaceName[s],space.name)==0)
      space_i = s;

  if(space_i == idfx::prof.numSpaces) {
    strncpy(idfx::prof.spaceName[space_i],space.name,64);
    idfx::prof.numSpaces++;
  }
  idfx::prof.spaceSize[space_i] -= size;

}
