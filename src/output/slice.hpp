// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/output/slice.hpp
//
// Last modified : 05/2024
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2023 - 2024)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************
// Other contributors :
// - Gaylor Wafflard <gaylor.wafflard@univ-grenoble-alpes.fr> (IPAG / CNRS - 2024)
// ***********************************************************************************

#ifndef OUTPUT_SLICE_HPP_
#define OUTPUT_SLICE_HPP_

#include <memory>
#include <map>
#include <string>
#include "idefix.hpp"
#include "dataBlock.hpp"
#include "input.hpp"

class Grid;
class SubGrid;
class Vtk;

using UserDefVariablesFunc = void (*) (DataBlock &,
                                      std::map<std::string,IdefixHostArray3D<real>> &);

class Slice {
 public:
  Slice(Input &, DataBlock &, int, SliceType, int, real, real);
  void CheckForWrite(DataBlock &, bool = false);
  void EnrollUserDefVariables(std::map<std::string,IdefixHostArray3D<real>>);
  void EnrollUserDefFunc(UserDefVariablesFunc);
  real slicePeriod = 0.0;
  real sliceLast = 0.0;
 private:
  bool containsX0;
  SliceType type;
  int direction;
  std::unique_ptr<SubGrid> subgrid;
  std::unique_ptr<DataBlock> sliceData;
  std::unique_ptr<Vtk> vtk;
  bool haveUserDefinedVariables{false};
  std::map<std::string, IdefixHostArray3D<real>> variableMap;
  std::map<std::string,IdefixHostArray3D<real>> userDefVariableMap;
  UserDefVariablesFunc userDefVariablesFunc{NULL};
  #ifdef WITH_MPI
    MPI_Comm avgComm;  // Communicator for averages
  #endif
};

#endif // OUTPUT_SLICE_HPP_
