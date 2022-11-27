#include "idefix.hpp"
#include "setup.hpp"


// Default constructor


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
              real r = d.x[IDIR](i);

                // Sedov Blast Wave Following Stone+2018, 3.4.4
                d.Vc(RHO,k,j,i) = 1.0;
                d.Vc(VX1,k,j,i) = 0.0;
                d.Vc(VX2,k,j,i) = 0.0;
                d.Vc(VX3,k,j,i) = 0.0;

                d.Vc(PRS,k,j,i) = 0.01;

                if(r<0.01) {
                  d.Vc(PRS,k,j,i) = 1.6e5;
                }


            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}
