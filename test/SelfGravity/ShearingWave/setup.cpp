#include "idefix.hpp"
#include "setup.hpp"

#ifndef ISOTHERMAL
static real gammaIdeal;
#endif
static real omega;
static real shear;
static real Gnewt;

//#define STRATIFIED

void BodyForce(DataBlock &data, const real t, IdefixArray4D<real> &force) {
  idfx::pushRegion("BodyForce");
  IdefixArray1D<real> x = data.x[IDIR];
  IdefixArray1D<real> z = data.x[KDIR];

  // GPUS cannot capture static variables
  real omegaLocal=omega;
  real shearLocal =shear;

  idefix_for("BodyForce",
              data.beg[KDIR] , data.end[KDIR],
              data.beg[JDIR] , data.end[JDIR],
              data.beg[IDIR] , data.end[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                EXPAND(
                force(IDIR,k,j,i) = -2.0*omegaLocal*shearLocal*x(i);,
                force(JDIR,k,j,i) = ZERO_F;,

                #ifdef STRATIFIED
                  force(KDIR,k,j,i) = - omegaLocal*omegaLocal*z(k);
                #else
                  force(KDIR,k,j,i) = ZERO_F;
                #endif
                );
      });


  idfx::popRegion();
}


void ComputeUserVars(DataBlock &data, UserDefVariablesContainer &variables) {
  Kokkos::deep_copy(variables["phiP"], data.gravity->phiP);
}
// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  #ifndef ISOTHERMAL
  gammaIdeal=data.hydro->eos->GetGamma();
  #endif

  // Get rotation rate along vertical axis
  omega=input.Get<real>("Hydro","rotation",0);
  shear=input.Get<real>("Hydro","shearingBox",0);
  Gnewt=input.Get<real>("Gravity","gravCst",0);

  // Add our userstep to the timeintegrator
  data.gravity->EnrollBodyForce(BodyForce);

  //output.EnrollUserDefVariables(&ComputeUserVars);

}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);
    real x,y,z = 0;

    // setup from Paardekooper 2012
    real kx = -4.*M_PI;
    real ky = 2.*M_PI;
    real sigma_init = 0.0005;
    real Sigma0 = 1./40.;
    real Q = 1.;

    real kappa = omega;

    real cs = M_PI*Sigma0*Q/kappa*Gnewt; // G = 1
    std::cout<< "cs = " << cs << std::endl;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                EXPAND(
                x=d.x[IDIR](i);,
                y=d.x[JDIR](j);,
                z=d.x[KDIR](k);
                );
                d.Vc(RHO,k,j,i) = Sigma0*(1.+  sigma_init*cos(kx*x+ky*y));
#ifndef ISOTHERMAL
                d.Vc(PRS,k,j,i) = d.Vc(RHO,k,j,i)*cs*cs;
#endif
                EXPAND(
                d.Vc(VX1,k,j,i) = 0.0;, //sigma_init/Sigma0*cos(kx*x+ky*y);,
                d.Vc(VX2,k,j,i) = shear*x;,// + sigma_init/Sigma0*cos(kx*x+ky*y);,
                d.Vc(VX3,k,j,i) = 0.0;
                );

            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}
