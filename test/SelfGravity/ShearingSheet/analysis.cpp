#include "analysis.hpp"
#include "fluid.hpp"
#include "idefix.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

Analysis::Analysis() {
  }

/* **************************************************************** */

void Analysis::PerformAnalysis(DataBlock &data) {
  idfx::pushRegion("Analysis::PerformAnalysis");
  data.DumpToFile("test");
  idfx::popRegion();
}
