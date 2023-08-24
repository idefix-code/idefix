#include "idefix.hpp"
#include "setup.hpp"


// Default constructor

real amplitude;
void InternalBoundary(DataBlock& data, const real t) {
  IdefixArray4D<real> Vc = data.hydro->Vc;
  idefix_for("InternalBoundary",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                // Cancel any motion that could be happening
            Vc(VX1,k,j,i) = 0.0;
                Vc(VX2,k,j,i) = 0.0;
                Vc(VX3,k,j,i) = 0.0;
              });
}

void UserDefBoundary(Hydro *hydro, int dir, BoundarySide side, real t) {
  auto *data = hydro->data;
  real B0 = 0.01;
  IdefixArray4D<real> Vc = data->hydro->Vc;
  IdefixArray4D<real> Vs = data->hydro->Vs;
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  IdefixArray1D<real> x3 = data->x[KDIR];

  IdefixArray1D<real> x1l = data->xl[IDIR];
  IdefixArray1D<real> x2l = data->xl[JDIR];
  IdefixArray1D<real> x3l = data->xl[KDIR];
  if (dir==IDIR) {
    int ibeg, iend;
    if (side == left) {
      ibeg = 0;
      iend = data->beg[IDIR];
      idefix_for("InnerBoundary", 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], ibeg, iend,
             KOKKOS_LAMBDA (int k, int j, int i) {
                 Vc(RHO,k,j,i) = 1.;
                 Vc(PRS,k,j,i) = 1.;
                 Vc(VX1,k,j,i) = ZERO_F;

                 Vs(BX1s,k,j,i) = ZERO_F;
                 Vs(BX2s,k,j,i) = ZERO_F;
                 Vs(BX3s,k,j,i) = B0*SIN(x2(j));
      });
    }
    if (side == right) {
      ibeg = data->end[IDIR];
      iend = data->np_tot[IDIR];
      idefix_for("InnerBoundary", 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], ibeg, iend,
             KOKKOS_LAMBDA (int k, int j, int i) {
                 Vc(RHO,k,j,i) = 1.;
                 Vc(PRS,k,j,i) = 1.;
                 Vc(VX1,k,j,i) = ZERO_F;

                 Vs(BX1s,k,j,i) = ZERO_F;
                 Vs(BX2s,k,j,i) = ZERO_F;
                 Vs(BX3s,k,j,i) = B0*SIN(x2(j));
      });
    }
  }
}

// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  data.hydro->EnrollUserDefBoundary(&UserDefBoundary);
  data.hydro->EnrollInternalBoundary(&InternalBoundary);

  amplitude = input.Get<real>("Setup","amplitude",0);
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
  // Create a host copy
  DataBlockHost d(data);
  real B0 = 0.01;
  IdefixArray4D<real> Vc = data.hydro->Vc;
  IdefixArray4D<real> Vs = data.hydro->Vs;
  IdefixArray1D<real> x1 = data.x[IDIR];
  IdefixArray1D<real> x2 = data.x[JDIR];
  IdefixArray1D<real> x3 = data.x[KDIR];
  IdefixArray1D<real> x1l = data.xl[IDIR];
  IdefixArray1D<real> x2l = data.xl[JDIR];
  IdefixArray1D<real> x3l = data.xl[KDIR];

  for(int k = 0; k < d.np_tot[KDIR] ; k++) {
    for(int j = 0; j < d.np_tot[JDIR] ; j++) {
      for(int i = 0; i < d.np_tot[IDIR] ; i++) {

          //compute an eigen mode for the temperature as a product of r, theta and phi eigen modes with l=2, m=1
          d.Vc(RHO,k,j,i) = (3./pow(x1(i),2.) - 1.)*SIN(x1(i))/x1(i) - 3*COS(x1(i))/pow(x1(i),2.); //jl spherical function
          d.Vc(RHO,k,j,i) *= -3.*COS(x2(j))*SIN(x2(j)); //Plm(cos(theta))
          d.Vc(RHO,k,j,i) *= SIN(x3(k)); //sin(m*phi)
          d.Vc(RHO,k,j,i) *= amplitude;
          d.Vc(RHO,k,j,i) += 1.;
          d.Vc(RHO,k,j,i) = 1./d.Vc(RHO,k,j,i);
          d.Vc(VX1,k,j,i) = ZERO_F;
          d.Vc(VX2,k,j,i) = ZERO_F;
          d.Vc(VX3,k,j,i) = ZERO_F;
          d.Vc(PRS,k,j,i) = 1.0;

          d.Vs(BX1s,k,j,i) = ZERO_F;
          d.Vs(BX2s,k,j,i) = ZERO_F;
          d.Vs(BX3s,k,j,i) = B0*SIN(x2(j));
      }
    }
  }

  // Send it all, if needed
  d.SyncToDevice();
}
