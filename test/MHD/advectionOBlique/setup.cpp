#include "idefix.hpp"
#include "setup.hpp"

/*********************************************/
/**
Customized random number generator
Allow one to have consistant random numbers
generators on different architectures.
**/
/*********************************************/

const int advectionDirection        =   IDIR;
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
    IdefixHostArray4D<real> A = IdefixHostArray4D<real>("VectorPotential",3,d.np_tot[KDIR],d.np_tot[JDIR],d.np_tot[IDIR]);

    real x,y;
    real B0=1e-6;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                switch(advectionDirection) {
                    case(IDIR):
                        x=d.x[IDIR](i);
                        y=d.x[JDIR](j);
                        break;
                    case(JDIR):
                        x=d.x[JDIR](j);
                        y=d.x[KDIR](k);
                        break;
                    case(KDIR):
                        IDEFIX_ERROR("can't do that direction");
                        break;
                }

                d.Vc(RHO,k,j,i) = 1.0;
                d.Vc(PRS,k,j,i) = 1.0;
                d.Vc(VX1,k,j,i) = 0.0;
                d.Vc(VX2,k,j,i) = 0.0;
                d.Vc(VX3,k,j,i) = 0.0;
                d.Vs(BX1s,k,j,i) = 0.0;
                d.Vs(BX2s,k,j,i) = 0.0;
                d.Vs(BX3s,k,j,i) = 0.0;

                A(IDIR,k,j,i) = 0.0;
                A(JDIR,k,j,i) = 0.0;
                A(KDIR,k,j,i) = B0*sin(2.0*M_PI*(x+y))/(2.0*M_PI);

                d.Vc(VX1+advectionDirection,k,j,i) = 1.0;
                d.Vc(VX1+advectionDirection+1,k,j,i) = 1.0;

            }
        }
    }

    // Create the staggered field from the magnetic potential
    d.MakeVsFromAmag(A);

    // Send it all, if needed
    d.SyncToDevice();
}
