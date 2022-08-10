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
    real x,y;

    real r0 = 0.1;
    real r1 = 0.115;
    real omega = 20;
    real B0 = 5.0/sqrt(4.0*M_PI);
    real vphi,f,r;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                x=d.x[IDIR](i);
                y=d.x[JDIR](j);

                r=pow(x*x+y*y,0.5);
                f=(r1-r)/(r1-r0);

                if(r<r0) {
                    d.Vc(RHO,k,j,i) = 10.0;
                    vphi = omega*r;
                }
                else if(r<=r1) {
                    d.Vc(RHO,k,j,i) = 1.0+9.0*f;
                    vphi = f*omega*r0;
                }
                else {
                    d.Vc(RHO,k,j,i) = 1.0;
                    vphi = 0.0;
                }
                d.Vc(PRS,k,j,i) = ONE_F;
                d.Vc(VX1,k,j,i) = -y/r*vphi;
                d.Vc(VX2,k,j,i) =x/r*vphi;


                d.Vs(BX1s,k,j,i) = B0;;
                d.Vs(BX2s,k,j,i) = 0.0;

            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}
