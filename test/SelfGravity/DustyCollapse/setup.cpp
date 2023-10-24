#include "idefix.hpp"
#include "setup.hpp"
#include "npy.hpp"
#include <iomanip>

real MTOT; // total mass of spheric particle clump
real R0;   // size of spheric particle clump


void ComputeUserVars(DataBlock & data, UserDefVariablesContainer &variables) {
  // Make references to the user-defined arrays (variables is a container of IdefixHostArray3D)
  // Note that the labels should match the variable names in the input file
  IdefixHostArray3D<real> phiP = variables["phiP"];
  Kokkos::deep_copy(phiP, data.gravity->phiP);
}

Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  MTOT = input.Get<real>("Setup", "Mtot", 0);
  R0 = input.Get<real>("Setup", "R0", 0);

  output.EnrollUserDefVariables(&ComputeUserVars);
}


real rho_profile (real r) { // density profile of spheric particle clump
  real rho_part = MTOT*3.0/(4.0*M_PI*R0*R0*R0);
  return (r < R0 ? rho_part : rho_part/10000);
}



void Setup::InitFlow(DataBlock &data) {
  DataBlockHost d(data);
  d.SyncFromDevice();

  IdefixHostArray1D<real> r = d.x[IDIR];

  // init flow
  for(int k = 0; k < d.np_tot[KDIR] ; k++) {
    for(int j = 0; j < d.np_tot[JDIR] ; j++) {
      for(int i = 0; i < d.np_tot[IDIR] ; i++) {
        d.Vc(RHO,k,j,i) = 1e-6;
        EXPAND(\
        d.Vc(VX1,k,j,i) = ZERO_F; ,\
        d.Vc(VX2,k,j,i) = ZERO_F; ,\
        d.Vc(VX3,k,j,i) = ZERO_F; )

        d.dustVc[0](RHO,k,j,i) = rho_profile(r(i));
        EXPAND(\
        d.dustVc[0](VX1,k,j,i) = ZERO_F; ,\
        d.dustVc[0](VX2,k,j,i) = ZERO_F; ,\
        d.dustVc[0](VX3,k,j,i) = ZERO_F; )

      }
    }
  }



  d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}
