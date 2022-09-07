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
    const int a    =    16807;
    const int m =    2147483647;
    static int in0 = 13763 + 2417*idfx::prank;
    int q;

    /* find random number  */
    q= (int) fmod((double) a * in0, m);
    in0=q;

    return((real) ((double) q/(double)m));
}


// User-defined boundaries
void UserdefBoundary(DataBlock& data, int dir, BoundarySide side, real t) {
    if( (dir==IDIR) && (side == left)) {
        IdefixArray4D<real> Vc = data.hydro.Vc;
        IdefixArray4D<real> Vs = data.hydro.Vs;
        IdefixArray1D<real> x1 = data.x[IDIR];

        int ighost = data.beg[IDIR];

        data.hydro.boundary.BoundaryFor("UserDefBoundary", dir, side,
            KOKKOS_LAMBDA (int k, int j, int i) {
                Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost);
                Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost);
                Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost) * sqrt(x1(i)/x1(ighost));
                Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost) * sqrt(x1(i)/x1(ighost));
                Vc(VX3,k,j,i) = Vc(VX3,k,j,ighost) * sqrt(x1(i)/x1(ighost));
            });
        data.hydro.boundary.BoundaryForX2s("UserDefBoundaryBX2s", dir, side,
            KOKKOS_LAMBDA (int k, int j, int i) {
                Vs(BX2s,k,j,i) = Vs(BX2s,k,j,ighost);
        });
        data.hydro.boundary.BoundaryForX3s("UserDefBoundaryBX3s", dir, side,
            KOKKOS_LAMBDA (int k, int j, int i) {
                Vs(BX3s,k,j,i) = Vs(BX3s,k,j,ighost);
            });
    }
}

void Potential(DataBlock& data, const real t, IdefixArray1D<real>& x1, IdefixArray1D<real>& x2, IdefixArray1D<real>& x3, IdefixArray3D<real>& phi) {

    idefix_for("Potential",0,data.np_tot[KDIR], 0, data.np_tot[JDIR], 0, data.np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
        phi(k,j,i) = -1.0/x1(i);
    });

}

// Default constructor


// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
    // Set the function for userdefboundary
    data.hydro.EnrollUserDefBoundary(&UserdefBoundary);
    data.gravity.EnrollPotential(&Potential);
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
                const real r=d.x[IDIR](i);

                d.Vc(RHO,k,j,i) = 1.0;
                d.Vc(PRS,k,j,i) = 1.0e-2;
                d.Vc(VX1,k,j,i) = 0.0;
                d.Vc(VX2,k,j,i) = pow(r,-0.5);
                d.Vc(VX3,k,j,i) = 1e-2*(0.5-randm());

                d.Vs(BX1s,k,j,i) = 0.0;
                d.Vs(BX2s,k,j,i) = 0.0;
                d.Vs(BX3s,k,j,i) = 1e-2*d.Vc(VX2,k,j,i);
            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}




// Do a specifically designed user step in the middle of the integration
void ComputeUserStep(DataBlock &data, real t, real dt) {

}
