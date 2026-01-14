#include "idefix.hpp"
#include "setup.hpp"


// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
      for(int j = 0; j < d.np_tot[JDIR] ; j++) {
        for(int i = 0; i < d.np_tot[IDIR] ; i++) {
          real x = d.x[IDIR](i);
          real y = d.x[JDIR](j);
          real z = d.x[KDIR](k);

          // Uniform density, no velocity
          d.Vc(RHO,k,j,i) = 1.0;
          d.Vc(VX1,k,j,i) = 0.0;
          d.Vc(VX2,k,j,i) = 0.0;
          d.Vc(VX3,k,j,i) = 0.0;

          #ifndef ISOTHERMAL
          d.Vc(PRS,k,j,i) = 1.0;
          #endif
        }
      }
    }
    // Send it all, if needed
    d.SyncToDevice();
}

