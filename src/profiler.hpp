// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef PROFILER_HPP_
#define PROFILER_HPP_

#include <map>
#include <mutex>  // NOLINT [build/c++11]
#include <string>

namespace idfx {

struct SpaceHandle {
  char name[64];
};

// Region is a helper class to Profiler
// it used to generate a tree of the regions encountered while running
// and produce a performance report.
class Region {
 public:
  Region(Region *parent, std::string name, int level);
  Region();
  ~Region();
  void Start();
  void Stop();
  void Show(double );
  Region* GetChild(std::string name);
  double GetTimer();
  static bool Compare(Region *, Region *);
  bool isLeaf{true};
  std::string name;
  std::map<std::string, Region*> children;
  Region *parent;
  int level;
 private:
  Kokkos::Timer timer;
  double myTime{0};
  int64_t nCalls{0};
};


class Profiler {
 public:
  void Init();
  void Show();
  void EnablePerformanceProfiling();
  int numSpaces;
  int64_t spaceSize[16];
  int64_t spaceMax[16];
  char spaceName[16][64];
  std::mutex m;

  bool perfEnabled{false};
  Region rootRegion;
  Region *currentRegion;
};



}// namespace idfx

#endif // PROFILER_HPP_
