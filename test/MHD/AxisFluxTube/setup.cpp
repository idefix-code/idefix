#include "idefix.hpp"
#include "setup.hpp"



real Rtorus;
real Ztorus;
real Rin;


// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  Rtorus = input.GetReal("Setup","Rtorus",0);
  Ztorus = input.GetReal("Setup","Ztorus",0);
  Rin = input.GetReal("Setup","Rin",0);
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
                real x1=d.x[IDIR](i);
                real x2=d.x[JDIR](j);
                real x3=d.x[KDIR](k);

                real R = x1*sin(x2);      /* Cylindrical radius */
                real Z = x1*cos(x2);      /* Vertical cylindrical coordinate */

                // Vector components in cartesian coordinates

                real ex_r=cos(x3)*sin(x2);
                real ex_t=cos(x3)*cos(x2);
                real ex_p=-sin(x3);

                real ey_r=sin(x3)*sin(x2);
                real ey_t=sin(x3)*cos(x2);
                real ey_p=cos(x3);

                real ez_r=cos(x2);
                real ez_t=-sin(x2);
                real ez_p=0.0;


                d.Vc(RHO,k,j,i) = 1.0;
                d.Vc(PRS,k,j,i) = 1.0;
                d.Vc(VX1,k,j,i) = (ex_r+ey_r)/sqrt(2);
                d.Vc(VX2,k,j,i) = (ex_t+ey_t)/sqrt(2);
                d.Vc(VX3,k,j,i) = (ex_p+ey_p)/sqrt(2);

                real bphi = 1.0 - (pow(R-Rtorus,2.0) + pow(Z-Ztorus,2.0)) / Rin;
                if(bphi<0.0) bphi = 0.0;
                d.Vs(BX1s,k,j,i) = 0.0;
                d.Vs(BX2s,k,j,i) = 0.0;
                d.Vs(BX3s,k,j,i) = 1.0e-8*bphi;


            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}
