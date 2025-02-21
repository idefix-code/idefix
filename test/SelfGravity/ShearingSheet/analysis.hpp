#ifndef ANALYSIS_HPP_
#define ANALYSIS_HPP_

#include "dataBlock.hpp"
#include "dataBlockHost.hpp"
#include "grid.hpp"
#include "idefix.hpp"
#include "input.hpp"
#include "output.hpp"
#include <fstream>
#include <iostream>

class Analysis {
public:
  // Constructor from Setup arguments
  Analysis();
  void PerformAnalysis(DataBlock &);
};

#endif // ANALYSIS_HPP__
