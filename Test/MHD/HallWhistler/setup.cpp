#include "idefix.hpp"
#include "setup.hpp"

/*********************************************/
/**
Customized random number generator
Allow one to have consistant random numbers
generators on different architectures.
**/
/*********************************************/

int mode;

// Default constructor
Setup::Setup() {}

// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Hydro &hydro) {
	mode = input.GetInt("Setup","mode",0);
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);
    real x,y,z;
		real kx,ky,kz;
		real lambdax,lambday,lambdaz;
		real B0x, B0y, B0z;
		real qx,qy,qz;

		// Setup an inclined field, not-grid alined, of norm 1
		/* 3D inclined */
		B0x=1/sqrt(14);
		B0y=2*B0x;
		B0z=3*B0x;

		/* 2D inclined */
		/*
		B0x=1/sqrt(5);
		B0z=2/sqrt(5);
		B0y=0.0;
		*/
		
		// inclined wavevector
		kx = 2.0*M_PI*B0x*mode;
		ky = 2.0*M_PI*B0y*mode;
		kz = 2.0*M_PI*B0z*mode;

		// Perturbation, perpendicular to the Field
		qx=2/sqrt(5);
		qz=-1/sqrt(5);
		qy=0.0;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                x=d.x[IDIR](i);
                y=d.x[JDIR](j);
                z=d.x[KDIR](k);

                d.Vc(RHO,k,j,i) = 1.0;
#if HAVE_ENERGY
                d.Vc(PRS,k,j,i) = 1.0;
#endif
                d.Vc(VX1,k,j,i) = qx*sin(kx*x+ky*y+kz*z);
								d.Vc(VX2,k,j,i) = qy*sin(kx*x+ky*y+kz*z);
								d.Vc(VX3,k,j,i) = qz*sin(kx*x+ky*y+kz*z);

								d.Vs(BX1s,k,j,i) = B0x;
#if DIMENSIONS >= 2
								d.Vs(BX2s,k,j,i) = B0y;
#else
								d.Vc(BX2,k,j,i) = B0y;
#endif
#if DIMENSIONS == 3
 								d.Vs(BX3s,k,j,i) = B0z;
#else
								d.Vc(BX3,k,j,i) = B0z;
#endif

            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void Setup::MakeAnalysis(DataBlock & data, real t) {

}



// Do a specifically designed user step in the middle of the integration
void ComputeUserStep(DataBlock &data, real t, real dt) {

}
