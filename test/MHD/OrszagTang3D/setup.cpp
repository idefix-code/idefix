#include "idefix.hpp"
#include "setup.hpp"

/*********************************************/
/**
Customized random number generator
Allow one to have consistant random numbers
generators on different architectures.
**/
/*********************************************/

Output* myOutput;
int outnum;
// Analysis function

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
    real x,y,z;
    IdefixHostArray4D<real> Ve;

    #ifndef EVOLVE_VECTOR_POTENTIAL
    Ve = IdefixHostArray4D<real>("Potential vector",3, d.np_tot[KDIR]+1, d.np_tot[JDIR]+1, d.np_tot[IDIR]+1);
    #else
    Ve = d.Ve;
    #endif

    bool haveTracer = data.hydro->haveTracer;

    real B0=1.0/sqrt(4.0*M_PI);

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                x=d.x[IDIR](i);
                y=d.x[JDIR](j);
                z=d.x[KDIR](k);

                d.Vc(RHO,k,j,i) = 25.0/(36.0*M_PI);
                d.Vc(PRS,k,j,i) = 5.0/(12.0*M_PI);
                d.Vc(VX1,k,j,i) = -sin(2.0*M_PI*y);
                d.Vc(VX2,k,j,i) = sin(2.0*M_PI*x)+cos(2.0*M_PI*z);
                d.Vc(VX3,k,j,i) = cos(2.0*M_PI*x);

                real xl=d.xl[IDIR](i);
                real yl=d.xl[JDIR](j);
                real zl=d.xl[KDIR](k);
                Ve(IDIR,k,j,i) = B0/(2.0*M_PI)*(cos(2.0*M_PI*yl));
                Ve(JDIR,k,j,i) = B0/(2.0*M_PI)*sin(2.0*M_PI*xl);
                Ve(KDIR,k,j,i) = B0/(2.0*M_PI)*(
                                    cos(2.0*M_PI*yl) + cos(4.0*M_PI*xl)/2.0);

                if(haveTracer) {
                  d.Vc(TRG  ,k,j,i) = x>0.5?  1.0:0.0;
                  d.Vc(TRG+1,k,j,i) = z>0.5?  1.0:0.0;
                }
            }
        }
    }

    #ifndef EVOLVE_VECTOR_POTENTIAL
    d.MakeVsFromAmag(Ve);
    #endif
    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}



// Do a specifically designed user step in the middle of the integration
void ComputeUserStep(DataBlock &data, real t, real dt) {

}
