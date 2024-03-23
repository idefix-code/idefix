#include "idefix.hpp"
#include "setup.hpp"


static int componentDir;     // Which field component varies
static int spatialDir;  // along Which direction it varies
static int coarseningDir; // coarsening direction (only one in this test!)
static real advSpeed;

void CoarsenFunction(DataBlock &data) {
  IdefixArray2D<int> coarseningLevel = data.coarseningLevel[coarseningDir];
  int dir = 0;
  if(coarseningDir == componentDir) dir = spatialDir;
  else if(coarseningDir == spatialDir) dir = componentDir;
  IdefixArray1D<real> x = data.x[dir];

  // Choose which index we should use for x, knowing that coarseningLevel lacks the direction coarseningDir
  bool useSecond = (dir == KDIR);
  if(dir == JDIR && coarseningDir == KDIR) useSecond = true;
  idefix_for("set_coarsening", 0, coarseningLevel.extent(0), 0, coarseningLevel.extent(1),
      KOKKOS_LAMBDA(int j,int i) {
        int idx = i;
        if(useSecond) idx = j;
        coarseningLevel(j,i) = 1 + static_cast<int>(6*(0.5-fabs(x(idx))));
        });
}


// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  componentDir = input.Get<int>("Setup", "direction", 0);
  spatialDir = input.Get<int>("Setup", "direction", 1);
  if(componentDir == spatialDir) {
    IDEFIX_ERROR("The field component and the direction dependence cannot be identical");
  }
  advSpeed = input.GetOrSet<real>("Setup","advectionSpeed",0,0.0);
  if(data.haveGridCoarsening) {
    data.EnrollGridCoarseningLevels(&CoarsenFunction);
    coarseningDir = -1;
    int i = 0;
    // detect coarsening direction
    while(coarseningDir < 0) {
      if(data.coarseningDirection[i]) coarseningDir = i;
      i++;
    }
  }

}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);

    int idx=0;
    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
      if(spatialDir == KDIR) idx = k;
      for(int j = 0; j < d.np_tot[JDIR] ; j++) {
        if(spatialDir == JDIR) idx = j;
        for(int i = 0; i < d.np_tot[IDIR] ; i++) {
          if(spatialDir == IDIR) idx = i;

          real x=d.x[spatialDir](idx);

          d.Vc(RHO,k,j,i) = 1.0;
          d.Vc(VX1,k,j,i) = 0.0;
          d.Vc(VX2,k,j,i) = 0.0;
          d.Vc(VX3,k,j,i) = 0.0;

          d.Vs(BX1s,k,j,i) = 0.0;
          d.Vs(BX2s,k,j,i) = 0.0;
          d.Vs(BX3s,k,j,i) = 0.0;

          d.Vs(componentDir,k,j,i) = fabs(x) < 0.1 ? 0.1 : 0.0;

          real Pmag = 0;
          for(int n = 0 ; n < DIMENSIONS ; n++) {
            Pmag += 0.5 * d.Vs(n,k,j,i) * d.Vs(n,k,j,i);
          }

          d.Vc(PRS,k,j,i) = 1.0 - Pmag;

          d.Vc(VX1+spatialDir,k,j,i) = advSpeed;
            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}
