// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/utils/dumpImage.hpp
//
// Last modified : 05/2024
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2021 - 2024)
// - Gaylor Wafflard <gaylor.wafflard@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2024)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef UTILS_DUMPIMAGE_HPP_
#define UTILS_DUMPIMAGE_HPP_
#include <string>
#include <map>
#include <set>
#include <vector>
#include "idefix.hpp"
#include "dump.hpp"

class DataBlock;

class DumpImage {
 public:
  DumpImage(std::string, DataBlock *, bool = false); // By default, no domain decomposition

  // Return a replicated array stored in the dump. The full (global) array is
  // returned on every rank, regardless of the domain decomposition.
  // Example: auto arr = image.GetSerialArray<double, 1>("myArray");
  template<typename T, int Rank>
  auto GetSerialArray(std::string const& name);

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

  // Replicated 1D/2D fields retained from the file, as flat byte buffers.
  // Populated during the read so that GetSerialArray() can serve them later.
  std::map<std::string, std::vector<char>> serialArrays;
  std::map<std::string, std::vector<int>> serialArrayShape;
  std::map<std::string, dump_io::DataType> serialArrayType;
};

template<typename T, int Rank>
auto DumpImage::GetSerialArray(std::string const& name) {
  if (serialArrays.count(name) == 0) {
    IDEFIX_ERROR("DumpImage: array " + name + " was not found in the dump");
  }
  std::vector<int> const& shape = serialArrayShape.at(name);
  if (static_cast<int>(shape.size()) != Rank) {
    IDEFIX_ERROR("DumpImage: array " + name + " does not have the requested rank");
  }
  if (serialArrayType.at(name) != dump_io::get_data_type<T>()) {
    IDEFIX_ERROR("DumpImage: array " + name + " does not have the requested type");
  }
  std::vector<char> const& buf = serialArrays.at(name);
  T const* src = reinterpret_cast<T const*>(buf.data());

  if constexpr (Rank == 1) {
    IdefixHostArray1D<T> out("DumpImageSerial", shape[0]);
    for (int i = 0 ; i < shape[0] ; i++) out(i) = src[i];
    return out;
  } else {
    IdefixHostArray2D<T> out("DumpImageSerial", shape[0], shape[1]);
    int n = 0;
    for (int i = 0 ; i < shape[0] ; i++) {
      for (int j = 0 ; j < shape[1] ; j++) out(i, j) = src[n++];
    }
    return out;
  }
}

#endif // UTILS_DUMPIMAGE_HPP_
