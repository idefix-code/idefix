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
    // Locally defined vector potential
    IdefixHostArray4D<real> A = IdefixHostArray4D<real>("Setup_VectorPotential", 3, data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);

    // Create a host copy
    DataBlockHost d(data);
    [[maybe_unused]] real x,y,z,r1,r2,r;

    real A0=1.0e-3;
    real R0=0.3;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {

                x=d.x[IDIR](i);
                y=d.x[JDIR](j);
                z=d.x[KDIR](k);

                r1 = sqrt(x*x + (z-1)*(z-1));
                r2 = sqrt(x*x + (z+1)*(z+1));
                r=sqrt(x*x + z*z);

                d.Vc(RHO,k,j,i) = 1.0;
                d.Vc(PRS,k,j,i) = 1.0;

                d.Vc(VX1,k,j,i) = 1.0;
                d.Vc(VX2,k,j,i) = 0.0;
                d.Vc(VX3,k,j,i) = 1.0;

                A(IDIR,k,j,i) = 0.0;
                //A(JDIR,k,j,i) = A0*(R0-r1)*(r1 <= R0);
                //A(JDIR,k,j,i) += A0*(R0-r2)*(r2 <= R0);
                A(JDIR,k,j,i) = A0*(R0-r)*(r <= R0);

                A(KDIR,k,j,i) = 0.0;

            }
        }
    }

    // Make the field from the vector potential
    d.MakeVsFromAmag(A);

    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}
