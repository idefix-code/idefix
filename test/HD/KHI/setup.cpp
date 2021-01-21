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


                // KHI
                d.Vc(RHO,k,j,i) = ONE_F;
                real yInterface = HALF_F*(1+0.05*(sin(0.5*M_PI*d.x[IDIR](i))+cos(4.0*M_PI*d.x[IDIR](i))));
                EXPAND(\
                d.Vc(VX1,k,j,i) = (d.x[JDIR](j) > yInterface) ? ONE_F : -ONE_F; ,\
                d.Vc(VX2,k,j,i) = ZERO_F; ,\
                d.Vc(VX3,k,j,i) = ZERO_F; )
#if HAVE_ENERGY
                d.Vc(PRS,k,j,i) = ONE_F;
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
