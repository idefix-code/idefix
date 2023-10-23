#include "idefix.hpp"
#include "setup.hpp"

// User-defined boundaries
void UserdefBoundary(Hydro *hydro, int dir, BoundarySide side, const real t) {
    auto *data = hydro->data;
    if(dir==IDIR) {
        IdefixArray4D<real> Vc = hydro->Vc;
        IdefixArray1D<real> th = data->x[JDIR];
        if(side==right) {
          int ighost = data->end[IDIR]-1;
          idefix_for("UserDefBoundaryX1End",
            0, data->np_tot[KDIR],
            0, data->np_tot[JDIR],
            data->end[IDIR], data->np_tot[IDIR],
                      KOKKOS_LAMBDA (int k, int j, int i) {
                          if((th(j)>M_PI/2.0) && th(j)<3.0*M_PI/2) {
                            // Incoming flow
                            Vc(RHO,k,j,i) = ONE_F;
                            Vc(VX1,k,j,i) = COS(th(j));
                            Vc(VX2,k,j,i) = -SIN(th(j));
                          } else {
                            // outcoming flow
                            Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost);
                            Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost);
                            if(Vc(VX1,k,j,ighost) > ZERO_F) {
                              Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost);
                            }else{
                              Vc(VX1,k,j,i) = ZERO_F;
                            }
                          }
                        });
      } else if(side==left) {
        int ighost = data->beg[IDIR];
        idefix_for("UserDefBoundaryX1Beg",
          0, data->np_tot[KDIR],
          0, data->np_tot[JDIR],
          0, ighost,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                      // Central cylinder
                      Vc(RHO,k,j,i) = ONE_F;
                      Vc(VX1,k,j,i) = -Vc(VX1,k,j,2*ighost-i-1);
                      Vc(VX2,k,j,i) = -Vc(VX2,k,j,2*ighost-i-1);;
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


    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                // Flow past a cylinder

                d.Vc(RHO,k,j,i) = ONE_F;
                d.Vc(VX1,k,j,i) = COS(d.x[JDIR](j));
                d.Vc(VX2,k,j,i) = -SIN(d.x[JDIR](j));

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
void MakeAnalysis(DataBlock & data) {

}
