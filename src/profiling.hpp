// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef PROFILING_HPP_
#define PROFILING_HPP_

#include <mutex>

namespace idfx {

struct SpaceHandle {
  char name[64];
};

class Profiler {
 public:
  void Init();
  void Show();
  int numSpaces;
  int64_t spaceSize[16];
  int64_t spaceMax[16];
  char spaceName[16][64];
  std::mutex m;
  };

}

#endif // PROFILING_HPP_
