// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef PROFILER_HPP_
#define PROFILER_HPP_

#include <mutex>  // NOLINT [build/c++11]

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

}// namespace idfx

#endif // PROFILER_HPP_
