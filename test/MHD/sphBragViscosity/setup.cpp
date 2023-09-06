#include "idefix.hpp"
#include "setup.hpp"


real amplitude;
void UserDefBoundary(Hydro *hydro, int dir, BoundarySide side, real t) {
  auto *data = hydro->data;
  IdefixArray4D<real> Vc = data->hydro->Vc;
  IdefixArray4D<real> Vs = data->hydro->Vs;
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  IdefixArray1D<real> x3 = data->x[KDIR];
  IdefixArray1D<real> x1l = data->xl[IDIR];
  if (dir==IDIR) {
    int ibeg, iend, iref;
    if (side == left) {
        ibeg = 0;
        iend = data->beg[IDIR];
        iref = data->beg[IDIR];
    }
    if (side == right) {
        ibeg = data->end[IDIR];
        iend = data->np_tot[IDIR];
        iref = data->end[IDIR] - 1;
    }
    idefix_for("InnerBoundary", 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], ibeg, iend,
          KOKKOS_LAMBDA (int k, int j, int i) {
          Vc(RHO,k,j,i) = 1.;

          Vc(VX1,k,j,i) = -Vc(VX1,k,j,iref);
          Vc(VX2,k,j,i) = ZERO_F;
          Vc(VX3,k,j,i) = ZERO_F;
    });
  }
}

// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  // Set the function for userdefboundary
  data.hydro->EnrollUserDefBoundary(&UserDefBoundary);

  amplitude = input.Get<real>("Setup","amplitude",0);
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
  // Create a host copy
  DataBlockHost d(data);
  real B0 = 1e-5;

  for(int k = 0; k < d.np_tot[KDIR] ; k++) {
    real x3 = d.x[KDIR](k);
    real x3l = d.xl[KDIR](k);
    for(int j = 0; j < d.np_tot[JDIR] ; j++) {
      real x2 = d.x[JDIR](j);
      real x2l = d.xl[JDIR](j);
      for(int i = 0; i < d.np_tot[IDIR] ; i++) {
        real x1=d.x[IDIR](i);
        real x1l=d.xl[IDIR](i);

        d.Vc(RHO,k,j,i) = 1.0;
        d.Vc(VX1,k,j,i) = amplitude*std::cyl_bessel_j(1, x1)/x1;
//        d.Vc(VX1,k,j,i) = amplitude;
        d.Vc(VX2,k,j,i) = 0.0;
        d.Vc(VX3,k,j,i) = 0.0;

        d.Vs(BX1s,k,j,i) = B0/pow(x1l,2.);
        d.Vs(BX2s,k,j,i) = 0.0;
        d.Vs(BX3s,k,j,i) = 0.0;
      }
    }
  }
  // Send it all, if needed
  d.SyncToDevice();
}
