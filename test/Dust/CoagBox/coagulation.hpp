#include "idefix.hpp"
#include "dataBlock/dataBlock.hpp"
#include "fluid/fluid.hpp"

class Coagulation {
 public:
  // Constructor: allocate and load all coagulation tables
  Coagulation(DataBlock &data, int kernel, real eps, real coagCFL);

  IdefixArray1D<real> binWidth;
  IdefixArray3D<real> coagTabFlux;

  int kernel;
  real eps;
  real coagCFL;
};

extern Coagulation *CoagulationGlob;

void CoagSourceTerm(DataBlock &data, const real t, const real dtin);