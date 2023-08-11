#include "idefix.hpp"
#include "setup.hpp"


// User-defined boundaries
void UserdefBoundary(Hydro *hydro, int dir, BoundarySide side, real t) {
    auto *data = hydro->data;
    if( (dir==IDIR) && (side == left)) {
        IdefixArray4D<real> Vc = hydro->Vc;
        IdefixArray4D<real> Vs = hydro->Vs;
        IdefixArray1D<real> x1 = data->x[IDIR];

        int ighost = data->nghost[IDIR];
        idefix_for("UserDefBoundary",
          0, data->np_tot[KDIR],
          0, data->np_tot[JDIR],
          0, ighost,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost);
                        Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost);
                        if(Vc(VX1,k,j,ighost) < ZERO_F) Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost) * sqrt(x1(i)/x1(ighost));
                        else Vc(VX1,k,j,i) = ZERO_F;
                        Vc(VX2,k,j,i) = 1/ sqrt(x1(i));
                        Vc(VX3,k,j,i) = Vc(VX3,k,j,ighost) * sqrt(x1(i)/x1(ighost));
                        Vs(BX2s,k,j,i) = Vs(BX2s,k,j,ighost);
                        Vs(BX3s,k,j,i) = Vs(BX3s,k,j,ighost);

                    });
    }

}

void Hall(DataBlock& data, const real t, IdefixArray3D<real> &xH) {
  IdefixArray1D<real> x1 = data.x[IDIR];

  idefix_for("Hall",0,data.np_tot[KDIR], 0, data.np_tot[JDIR], 0, data.np_tot[IDIR],

      KOKKOS_LAMBDA (int k, int j, int i) {
        real f;
        if(x1(i)<1.5) f = 0.0;
        else if(x1(i) < 1.6) f = 10.0*(x1(i)-1.5);
        else f=1.0;

    real g=0.0;
    if(t>200.0 && t< 201.0) g=t-100.0;
          if(t>201.0) g=1.0;
        xH(k,j,i) = 0.1*g*pow(x1(i),-0.5)*f;
  });
}
// Default constructor


// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
    // Set the function for userdefboundary
    data.hydro->EnrollUserDefBoundary(&UserdefBoundary);
    data.hydro->EnrollHallDiffusivity(&Hall);
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);

    real beta=1e3;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                const real r = d.x[IDIR](i);

                d.Vc(RHO,k,j,i) = 1.0;
                d.Vc(PRS,k,j,i) = 1.0e-2;
                d.Vc(VX1,k,j,i) = 0.0;
                d.Vc(VX2,k,j,i) = pow(r,-0.5);
                d.Vc(VX3,k,j,i) = 1e-2*(0.5-idfx::randm());

                d.Vs(BX1s,k,j,i) = 0.0;
                d.Vs(BX2s,k,j,i) = 0.0;
                d.Vs(BX3s,k,j,i) = 1e-1*d.Vc(VX2,k,j,i)/sqrt(beta);
            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}
