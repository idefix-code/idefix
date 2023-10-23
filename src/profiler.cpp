// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <algorithm>
#include <iomanip>
#include <mutex>    // NOLINT [build/c++11]
#include <string>
#include <vector>

#include "idefix.hpp"
#include "profiler.hpp"

/////////////////////////////
// Kokkos Profiler hooks (needed to track memory allocation)
/////////////////////////////

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

///////////////////////////////////
// Profiler function definitions //
///////////////////////////////////
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

  // Init region


  idfx::popRegion();
}

void idfx::Profiler::Show() {
  double usedMemory{-1};

  // follow ISO/IEC 80000
  constexpr int nUnits = 5;
  const std::array<std::string,nUnits> units {"B", "KB", "MB", "GB", "TB"};

  for(int i=0; i < this->numSpaces ; i++) {
    usedMemory = this->spaceMax[i];
    int count{0};
    while(count < nUnits && usedMemory/1024 >= 1) {
      usedMemory /= 1024;
      ++count;
    }
    idfx::cout << "Profiler: maximum memory usage for " << this->spaceName[i];
    idfx::cout << " memory space: " << usedMemory << " " << units[count] << std::endl;
  }

  if(perfEnabled) {
    // Show performance results
    rootRegion.Stop();
    idfx::cout << "Profiler: performance results: " << std::endl;
    idfx::cout << "-------------------------------------------------------------------------------";
    idfx::cout << std::endl;
    idfx::cout << "<total time>  <% of total time>  <% of self time>  <number of calls>  <name>";
    idfx::cout << std::endl;
    idfx::cout << "-------------------------------------------------------------------------------";
    idfx::cout << std::endl;
    rootRegion.Show(rootRegion.GetTimer());
    idfx::cout << "-------------------------------------------------------------------------------";
    idfx::cout << std::endl;
    idfx::cout << "Profiler: end of performance profiling report." << std::endl;
  }
}

void idfx::Profiler::EnablePerformanceProfiling() {
  currentRegion = &rootRegion;
  rootRegion.Start();
  perfEnabled = true;
}


///////////////////////////////////
// Region functions definitions //
///////////////////////////////////


idfx::Region::Region(Region *parent, std::string name, int level) {
  this->parent = parent;
  this->name = name;
  this->level = level;
}

idfx::Region::Region() {
  this->parent = nullptr;
  this->name = std::string("Main");
  this->level = 0;
}

idfx::Region::~Region() {
  if(!isLeaf) {
    // Delete all the children manually, since these are pointers
    for( auto &it : children) {
      delete it.second;
    }
  }
}

void idfx::Region::Start() {
  this->nCalls++;
  this->timer.reset();
}

double idfx::Region::GetTimer() {
  return this->myTime;
}

void idfx::Region::Stop() {
  this->myTime += this->timer.seconds();
}

idfx::Region * idfx::Region::GetChild(std::string name) {
  isLeaf=false;
  if(this->children.find(name) == this->children.end()) {
    // No children found
    this->children[name] = new Region(this, name, this->level+1);
  }
  return this->children[name];
}

bool idfx::Region::Compare(Region * r1, Region * r2) {
  return r1->GetTimer() > r2->GetTimer();
}

void idfx::Region::Show(double totTime) {
  // Compute time of all the children
  double childTime = 0;

  if(!isLeaf) {
    for( auto &it : children) {
          childTime += it.second->GetTimer();
        }
  }
  for(int i = 0 ; i < (this->level) ; i++) {
    idfx::cout << "|   ";
  }

  idfx::cout << "|-> " << std::scientific << std::setprecision(2) << this->myTime << " sec  "
             << std::fixed << std::setprecision(1)
             << this->myTime/totTime*100 << "%  "
             << (this->myTime-childTime)/this->myTime*100 << "%  "
             << this->nCalls << "  "
             << this->name << std::endl;
  if(!isLeaf) {
    // Sort the children
    std::vector<Region*> sorted;
    for( auto &it : this->children) {
      sorted.push_back(it.second);
    }

    std::sort(sorted.begin(), sorted.end(), this->Compare);
    for( auto &it : sorted) {
          it->Show(totTime);
        }
  }
}
