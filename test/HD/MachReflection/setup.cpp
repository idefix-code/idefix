#include "idefix.hpp"
#include "setup.hpp"

// User-defined boundaries
void UserdefBoundary(Hydro* hydro, int dir, BoundarySide side, const real t) {
    const real alpha = 60./180.*M_PI;

    if( (dir==IDIR) && (side == left)) {
        IdefixArray4D<real> Vc = hydro->Vc;
        int ighost = hydro->data->nghost[IDIR];
        hydro->boundary->BoundaryFor("UserDefBoundaryX1Beg", dir, side,
                    KOKKOS_LAMBDA (int k, int j, int i) {

                        Vc(RHO,k,j,i) = 8.0;
                        Vc(VX1,k,j,i) =   8.25*sin(alpha);
                        Vc(VX2,k,j,i) = - 8.25*cos(alpha);
                        Vc(PRS,k,j,i) = 116.5;
                      });
    }

    if(dir==JDIR) {
        IdefixArray4D<real> Vc = hydro->Vc;
        IdefixArray1D<real> x1 = hydro->data->x[IDIR];
        IdefixArray1D<real> x2 = hydro->data->x[JDIR];
        if(side == left) {
            const int jend = hydro->data->beg[JDIR];
            hydro->boundary->BoundaryFor("UserDefBoundaryX2Beg", dir, side,
                        KOKKOS_LAMBDA (int k, int j, int i) {
                            if(x1(i) < 1.0/6.0) {
                              Vc(RHO,k,j,i) = 8.0;
                              Vc(VX1,k,j,i) =   8.25*sin(alpha);
                              Vc(VX2,k,j,i) = - 8.25*cos(alpha);
                              Vc(PRS,k,j,i) = 116.5;
                            } else {
                              Vc(RHO,k,j,i) =  Vc(RHO,k,2*jend-j-1,i);
                              Vc(VX1,k,j,i) =  Vc(VX1,k,2*jend-j-1,i);
                              Vc(VX2,k,j,i) = -Vc(VX2,k,2*jend-j-1,i);
                              Vc(PRS,k,j,i) =  Vc(PRS,k,2*jend-j-1,i);
                            }

                        });
            //return;
        } else if(side==right) {
            real xs = 10.0*t/sin(alpha) + 1.0/6.0 + 1.0/tan(alpha);

            hydro->boundary->BoundaryFor("UserDefBoundaryX2End", dir, side,
                        KOKKOS_LAMBDA (int k, int j, int i) {
                            if(x1(i) < xs) {
                              Vc(RHO,k,j,i) = 8.0;
                              Vc(VX1,k,j,i) =   8.25*sin(alpha);
                              Vc(VX2,k,j,i) = - 8.25*cos(alpha);
                              Vc(PRS,k,j,i) = 116.5;
                            } else {
                              Vc(RHO,k,j,i) =  1.4;
                              Vc(VX1,k,j,i) =  0.0;
                              Vc(VX2,k,j,i) =  0.0;
                              Vc(PRS,k,j,i) =  1.0;
                            }

                        });
        }
    }

}

// Default constructor


// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
    data.hydro->EnrollUserDefBoundary(&UserdefBoundary);
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);
    real alpha, xs, x1,x2;

    alpha = 1.0/3.0*M_PI;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                x1 = d.x[IDIR](i);
                x2 = d.x[JDIR](j);
                xs = 1.0/6.0 + x2/tan(alpha);
                if(x1>xs) {
                  d.Vc(RHO,k,j,i) = 1.4;
                  d.Vc(VX1,k,j,i) = 0.0;
                  d.Vc(VX2,k,j,i) = 0.0;
                  d.Vc(PRS,k,j,i) = 1.0;
                } else {
                  d.Vc(RHO,k,j,i) = 8.0;
                  d.Vc(VX1,k,j,i) = 8.25*sin(alpha);
                  d.Vc(VX2,k,j,i) = -8.25*cos(alpha);
                  d.Vc(PRS,k,j,i) = 116.5;
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



// Do a specifically designed user step in the middle of the integration
void ComputeUserStep(DataBlock &data, real t, real dt) {

}
