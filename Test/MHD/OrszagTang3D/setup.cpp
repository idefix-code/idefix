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
Setup::Setup() {}

// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, TimeIntegrator &tint) {

}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally 
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);
    real x,y,z;

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
                d.Vc(VX2,k,j,i) = sin(2.0*M_PI*x);
                d.Vc(VX3,k,j,i) = 0.0;
                
                d.Vs(BX1s,k,j,i) = -B0*sin(2.0*M_PI*y);
                d.Vs(BX2s,k,j,i) = B0*sin(4.0*M_PI*x);
                d.Vs(BX3s,k,j,i) = B0*(cos(2.0*M_PI*x)+sin(2.0*M_PI*y));

            }
        }
    }
    
    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void Setup::MakeAnalysis(DataBlock & data, real t) {

}

// User-defined boundaries
void Setup::SetUserdefBoundary(DataBlock& data, int dir, BoundarySide side, real t) {

}


// Do a specifically designed user step in the middle of the integration
void ComputeUserStep(DataBlock &data, real t, real dt) {

}
