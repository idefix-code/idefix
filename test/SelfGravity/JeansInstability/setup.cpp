#include <cmath>

#include "idefix.hpp"
#include "setup.hpp"

real lengthA, xc, sigma;
real rho0, gammaGlob, prs0, A;

// Default constructor


// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  lengthA = grid.xend[IDIR]-grid.xbeg[JDIR]; // Length of the active total domain
  xc = 0.5 * lengthA; // gaussian's center
  sigma = 0.1 * lengthA; // gaussian's variance

  // Parameters of the physical setup
  rho0 = 3.0; // Background density
  gammaGlob=data.hydro->eos->GetGamma(); // Input gamma
  prs0 = 1.0/gammaGlob; // Background pressure
  A = 0.0001; // Perturbation amplitude
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);

    // Grid parameters
    IdefixHostArray1D<real> x1 = d.x[IDIR]; // Position on the dataBlock

    // Mapping position wrt 0. and ensuring periodicity
    IdefixHostArray1D<real> y1 = IdefixHostArray1D<real> ("CenteredPosition", d.np_tot[IDIR]); // Position wrt to the new center 0.

    for(int i = 0; i < d.np_tot[IDIR] ; i++) {
      if((x1(i) - xc) > lengthA/2.) {
        y1(i) = x1(i) - xc - lengthA;
      }
      else if((x1(i) - xc) < -1.*lengthA/2.) {
        y1(i) = x1(i) - xc + lengthA;
      }
      else {
        y1(i) = x1(i) - xc;
      }
    }

    // Filling initial profiles
    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                d.Vc(RHO,k,j,i) = rho0 * (1. + A * exp(-1.*y1(i)*y1(i)/(2.*sigma*sigma)));
                d.Vc(VX1,k,j,i) = 0.;
                d.Vc(PRS,k,j,i) = prs0 * (1. + gammaGlob * A * exp(-1.*y1(i)*y1(i)/(2.*sigma*sigma)));
            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}
