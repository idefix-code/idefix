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

                d.Vc(RHO,k,j,i) = (d.x[IDIR](i)>HALF_F) ? 0.125 : 1.0;
                d.Vc(VX1,k,j,i) = ZERO_F;
#if HAVE_ENERGY
                d.Vc(PRS,k,j,i) = (d.x[IDIR](i)>HALF_F) ? 0.1 : 1.0;
#endif

            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}
