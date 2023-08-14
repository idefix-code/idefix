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
