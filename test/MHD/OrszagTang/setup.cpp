#include "idefix.hpp"
#include "setup.hpp"

/*********************************************/
/**
Customized random number generator
Allow one to have consistant random numbers
generators on different architectures.
**/
/*********************************************/


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

    bool haveTracer = data.hydro->haveTracer;

    real B0=1.0/sqrt(4.0*M_PI);

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                real x=d.x[IDIR](i);
                real y=d.x[JDIR](j);

                d.Vc(RHO,k,j,i) = 25.0/(36.0*M_PI);
                d.Vc(PRS,k,j,i) = 5.0/(12.0*M_PI);
                d.Vc(VX1,k,j,i) = -sin(2.0*M_PI*y);
                d.Vc(VX2,k,j,i) = sin(2.0*M_PI*x);
                #ifdef EVOLVE_VECTOR_POTENTIAL
                  x=d.xl[IDIR](i);
                  y=d.xl[JDIR](j);
                  d.Ve(AX3e,k,j,i) = B0/(2.0*M_PI)*(
                                      cos(2.0*M_PI*y) + cos(4.0*M_PI*x)/2.0);
                #else
                  d.Vs(BX1s,k,j,i) = -B0*sin(2.0*M_PI*y);
                  d.Vs(BX2s,k,j,i) = B0*sin(4.0*M_PI*x);
                #endif
                if(haveTracer) {
                  d.Vc(TRG,k,j,i) = x>0.5?  1.0:0.0;
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
