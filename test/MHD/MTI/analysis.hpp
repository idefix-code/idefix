// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file test/MHD/MTI/analysis.hpp
//
// Last modified : 10/2023
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG / CNRS - 2023)
// - and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef ANALYSIS_HPP_
#define ANALYSIS_HPP_

#include "idefix.hpp"
#include "input.hpp"
#include "output.hpp"
#include "grid.hpp"
#include "dataBlock.hpp"
#include "dataBlockHost.hpp"
#include <iostream>
#include <fstream>


class Analysis {
 public:
  // Constructor from Setup arguments
  Analysis(Input&, Grid&, DataBlock& , Output&, std::string);
  void ResetAnalysis();
  void PerformAnalysis(DataBlock &);

 private:
  double Average(const int, int[]);
  void WriteField(double);

  DataBlockHost *d;
  DataBlock *datain;
  Grid *grid;

  int precision;
  std::string filename;

  std::ofstream file;
};

#endif // ANALYSIS_HPP__
