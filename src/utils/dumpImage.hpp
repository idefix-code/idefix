// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef UTILS_DUMPIMAGE_HPP_
#define UTILS_DUMPIMAGE_HPP_
#include <string>
#include <map>
#include "idefix.hpp"
#include "dump.hpp"

class DataBlock;

class DumpImage {
 public:
  DumpImage(std::string, DataBlock *, bool = false); // By default, no domain decomposition

  int np_int[3];               // number of local points in each direction
  int np_glob[3];               // number of global points in each direction
  int geometry;                // geometry of the dump
  real time;                   // time at which the dump was created
  real centralMass;            // central mass when dump was created
  std::array<IdefixHostArray1D<real>,3> x;    // geometrical central points
  std::array<IdefixHostArray1D<real>,3> xr;   // cell right interface
  std::array<IdefixHostArray1D<real>,3> xl;   // cell left interface

  std::map<std::string,IdefixHostArray3D<real>> arrays;  // 3D arrays stored in the dump
 private:
  GridBox GetBox(DataBlock *);
};

#endif // UTILS_DUMPIMAGE_HPP_
