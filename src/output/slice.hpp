// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************


#ifndef OUTPUT_SLICE_HPP_
#define OUTPUT_SLICE_HPP_

#include <memory>
#include "idefix.hpp"
#include "dataBlock.hpp"
#include "input.hpp"

class Grid;
class SubGrid;
class Vtk;


class Slice {
 public:
  Slice(Input &, DataBlock &, int, SliceType, int, real, real);
  void CheckForWrite(DataBlock &);
  real slicePeriod = 0.0;
  real sliceLast = 0.0;
 private:
  bool containsX0;
  IdefixArray4D<real> Vc;
  SliceType type;
  int direction;
  std::unique_ptr<SubGrid> subgrid;
  std::unique_ptr<DataBlock> sliceData;
  std::unique_ptr<Vtk> vtk;
};

#endif // OUTPUT_SLICE_HPP_
