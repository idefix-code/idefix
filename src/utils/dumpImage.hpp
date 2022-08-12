// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef UTILS_DUMPIMAGE_HPP_
#define UTILS_DUMPIMAGE_HPP_
#include <string>
#include <map>
#include "idefix.hpp"
#include "dump.hpp"

class DumpImage {
 public:
  DumpImage(std::string, Output &);

  int np_int[3];               // number of points in each direction
  int geometry;                // geometry of the dump
  real time;                   // time at which the dump was created
  IdefixHostArray1D<real> x[3];    // geometrical central points
  IdefixHostArray1D<real> xr[3];   // cell right interface
  IdefixHostArray1D<real> xl[3];   // cell left interface

  std::map<std::string,IdefixHostArray3D<real>> arrays;  // 3D arrays stored in the dump
};

#endif // UTILS_DUMPIMAGE_HPP_
