#include <cmath>

#include "idefix.hpp"
#include "setup.hpp"

real delta, r0, rho0;

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
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);

    // Grid parameters
    IdefixHostArray1D<real> r = d.x[IDIR];
    IdefixHostArray1D<real> th = d.x[JDIR];
    IdefixHostArray1D<real> phi = d.x[KDIR];

    // Sphere's parameters
    r0 = 3.0; // Sphere's center position
    delta = 1.0; // Sphere's radius
    rho0 = 3.0/(4.*M_PI); // Sphere's density

    // Filling initial profiles
    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
              real x, y, z;
              x = r(i)*sin(th(j))*cos(phi(k));
              y = r(i)*sin(th(j))*sin(phi(k));
              z = r(i)*cos(th(j));
              real dr = sqrt((x-r0)*(x-r0) + y*y + z*z);
              d.Vc(RHO,k,j,i) = (dr<=delta) ? rho0 : ZERO_F;
            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}
