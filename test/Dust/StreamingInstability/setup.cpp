#include "idefix.hpp"
#include "setup.hpp"

static real omega;
static real shear;
real epsilon;
real chi;
real tauGlob;


#define  FILENAME    "timevol.dat"

//#define STRATIFIED
void PressureGradient(Hydro *hydro, const real t, const real dt) {
  auto Uc = hydro->Uc;
  auto Vc = hydro->Vc;
  DataBlock *data = hydro->data;
  real eps = epsilon;
  idefix_for("MySourceTerm",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                // Radial pressure gradient
                  Uc(MX1,k,j,i) += eps*Vc(RHO,k,j,i)*dt;
              });
}

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

                force(IDIR,k,j,i) = -2.0*omegaLocal*shearLocal*x(i);
                force(JDIR,k,j,i) = ZERO_F;
                force(KDIR,k,j,i) = ZERO_F;

      });


  idfx::popRegion();
}

// Analyse data to produce an output
void Analysis(DataBlock & data) {



  if(idfx::prank == 0) {
    std::ofstream f;
    f.open(FILENAME,std::ios::app);
    f.precision(10);
    int k=18;
    int j=18;
    int i=18;

    real vx = data.hydro->Vc(VX1,k,j,i);
    real vy = data.hydro->Vc(VX2,k,j,i);
    real ux = data.dust[0]->Vc(VX1,k,j,i);
    real uy = data.dust[0]->Vc(VX2,k,j,i);
    f << std::scientific << data.t << "\t" << vx << "\t" << vy << "\t"  << ux << "\t" << uy << std::endl;
    f.close();
  }

}

// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  // Get rotation rate along vertical axis
  omega=input.Get<real>("Hydro","rotation",0);
  shear=input.Get<real>("Hydro","shearingBox",0);

  tauGlob = input.Get<real>("Dust","drag",1);
  epsilon = input.Get<real>("Setup","epsilon",0);
  chi = input.Get<real>("Setup","chi",0);
  data.hydro->EnrollUserSourceTerm(&PressureGradient);
  // Add our userstep to the timeintegrator
  data.gravity->EnrollBodyForce(BodyForce);

}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);

    real taus = tauGlob*omega;
    real D = 1+chi;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                real x=d.x[IDIR](i);


                d.Vc(RHO,k,j,i) = 1.0;
                // Equations 7a-7b of Youdin & Johansen (2007) with eta = epsilon/(2Omega Vk)
                d.Vc(VX1,k,j,i) = chi * taus/((D*D+taus*taus))*epsilon/omega;
                d.Vc(VX2,k,j,i) = -(1+ chi*taus*taus/(D*D+taus*taus))*epsilon/(2*D*omega);
                d.Vc(VX1,k,j,i) += 1e-1*(idfx::randm()-0.5);
                d.Vc(VX2,k,j,i) += shear*x;
                d.Vc(VX3,k,j,i) = 0.0;

                d.dustVc[0](RHO,k,j,i) = chi;
                // Equations 7c-7d of Youdin & Johansen (2007) with eta = epsilon/(2Omega Vk)
                d.dustVc[0](VX1,k,j,i) = - taus/(D*D+taus*taus)*epsilon/omega;
                d.dustVc[0](VX2,k,j,i) = - (1-taus*taus/(D*D+taus*taus))*epsilon/(2*D*omega);

                d.dustVc[0](VX1,k,j,i) += 0*1e-1*(idfx::randm()-0.5);
                d.dustVc[0](VX2,k,j,i) += shear*x;
                d.dustVc[0](VX3,k,j,i) = 0.0;

            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}


// Analyse data to produce an output

void MakeAnalysis(DataBlock & data) {
}
