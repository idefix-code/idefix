#include <cmath>

#include "idefix.hpp"
#include "setup.hpp"

real inx0, iny0, inz0, inr0;

// Compute user variables which will be written in vtk files
void ComputeUserVars(DataBlock & data, UserDefVariablesContainer &variables) {
  // Force compute the gravity field
  idfx::cout << "ComputeUserVars: Computing total gravity field..." << std::endl;
  data.gravity->ComputeGravity(0);
  idfx::cout << "ComputeUserVars: Done." << std::endl;

  // Make references to the user-defined arrays (variables is a container of IdefixHostArray3D)
  // Note that the labels should match the variable names in the input file
  // Storing into variables
  Kokkos::deep_copy(variables["phiP"],data.gravity->phiP);
}


// Default constructor


// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  // Enroll output functions
  output.EnrollUserDefVariables(&ComputeUserVars);
  inx0 = input.Get<real>("Setup","x0",0);
  iny0 = input.Get<real>("Setup","y0",0);
  inz0 = input.Get<real>("Setup","z0",0);
  inr0 = input.Get<real>("Setup","r0",0);
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);

    real rho0 = 1.0; // Sphere's density
    real x0=inx0;
    real y0=iny0;
    real z0=inz0;
    real r0=inr0;
    // Filling initial profiles
    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
              real dx = d.x[IDIR](i)-x0;
              real dy = d.x[JDIR](j)-y0;
              real dz = d.x[KDIR](k)-z0;
              real dr = sqrt(dx*dx+dy*dy+dz*dz);
              d.Vc(RHO,k,j,i) = (dr<=r0) ? rho0 : ZERO_F;
            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}
