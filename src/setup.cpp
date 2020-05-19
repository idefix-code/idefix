#include "idefix.hpp"
#include "setup.hpp"

/*********************************************/
/**
Customized random number generator
Allow one to have consistant random numbers
generators on different architectures.
**/
/*********************************************/
real randm(void) {
	const int a	=	16807;
	const int m =	2147483647;
	static int in0 = 13763;
	int q;

	/* find random number  */
	q= (int) fmod((double) a * in0, m);
	in0=q;

	return((real) ((double) q/(double)m));
}


// Default constructor
Setup::Setup() {}

// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, TimeIntegrator &Tint) {

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
                
                /*
                // SHOCK TUBE
                d.Vc(RHO,k,j,i) = (d.x[IDIR](i)>HALF_F) ? 0.125 : 1.0;
                d.Vc(VX1,k,j,i) = ZERO_F;
#if HAVE_ENERGY 
                d.Vc(PRS,k,j,i) = (d.x[IDIR](i)>HALF_F) ? 0.1 : 1.0;
#endif
                */
                // KHI
                d.Vc(RHO,k,j,i) = ONE_F;
                EXPAND(\
                d.Vc(VX1,k,j,i) = (d.x[JDIR](j) > HALF_F) ? ONE_F : -ONE_F; ,\
                d.Vc(VX2,k,j,i) = randm()-HALF_F; ,\
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
void Setup::MakeAnalysis(DataBlock & data, real t) {

}

// User-defined boundaries
void Setup::SetUserdefBoundary(DataBlock& data, int dir, BoundarySide side, real t) {

}


// Do a specifically designed user step in the middle of the integration
void ComputeUserStep(DataBlock &data, real t, real dt) {

}