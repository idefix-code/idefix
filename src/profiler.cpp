// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <mutex>    // NOLINT [build/c++11]
#include "idefix.hpp"
#include "profiler.hpp"

// Kokkos Profiler hooks

extern "C" void kokkosp_allocate_data(const Kokkos_Profiling_SpaceHandle space,
  const char* label, const void* const ptr, const uint64_t size) {
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


extern "C" void kokkosp_deallocate_data(const Kokkos_Profiling_SpaceHandle space,
const char* label, const void* const ptr, const uint64_t size) {
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


void idfx::Profiler::Init() {
  idfx::pushRegion("Profiler::Init");

  this->numSpaces=0;
  for(int i=0; i < 16 ; i++) {
    this->spaceSize[i] = 0;
    this->spaceMax[i] = 0;
  }

  // enroll callback
  Kokkos::Tools::Experimental::set_allocate_data_callback(&kokkosp_allocate_data);
  Kokkos::Tools::Experimental::set_deallocate_data_callback(&kokkosp_deallocate_data);

  idfx::popRegion();
}

void idfx::Profiler::Show() {
  idfx::pushRegion("Profiler::Show");

  for(int i=0; i < this->numSpaces ; i++) {
    idfx::cout << "Profiler: maximum memory usage for " << this->spaceName[i] << " memory space: "
               <<  this->spaceMax[i]/(1024.0*1024.0) << " MB." << std::endl;
  }
  idfx::popRegion();
}
