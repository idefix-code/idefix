#include "idefix.hpp"
#include "setup.hpp"

static real gammaIdeal;
static real omega;
static real shear;

//#define STRATIFIED

// UserStep, here only gravity (vertical and radial)
void UserStep(Hydro *hydro, const real t, const real dt) {
    auto *data = hydro->data;
    Kokkos::Profiling::pushRegion("Setup::UserStep");
    IdefixArray4D<real> Uc = hydro->Uc;
    IdefixArray4D<real> Vc = hydro->Vc;
    IdefixArray1D<real> x = data->x[IDIR];
    IdefixArray1D<real> z = data->x[KDIR];

    // GPUS cannot capture static variables
    real omegaLocal=omega;
    real shearLocal =shear;

    idefix_for("UserSourceTerms",
      data->beg[KDIR], data->end[KDIR],
      data->beg[JDIR], data->end[JDIR],
      data->beg[IDIR], data->end[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
#ifdef STRATIFIED
            Uc(MX3,k,j,i) += - dt*omegaLocal*omegaLocal*z(k)*Vc(RHO,k,j,i);
#endif
            Uc(MX1,k,j,i) += - TWO_F*dt*omegaLocal*shearLocal*Vc(RHO,k,j,i)*x(i);
        });

    Kokkos::Profiling::popRegion();

}


// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
    gammaIdeal=data.hydro->eos->GetGamma();

    // Get rotation rate along vertical axis
    omega=input.Get<real>("Hydro","rotation",0);
    shear=input.Get<real>("Hydro","shearingBox",0);

    // Add our userstep to the timeintegrator
    data.hydro->EnrollUserSourceTerm(UserStep);
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);

    real B0z = 0.025;
    real B0y = 0.1;
    real cs2 = gammaIdeal*omega*omega;


    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                real x = d.x[IDIR](i);
                [[maybe_unused]] real z=d.x[KDIR](k);

#ifdef STRATIFIED
                d.Vc(RHO,k,j,i) = 1.0*exp(-z*z/(2.0));
#else
                d.Vc(RHO,k,j,i) = 1.0;
#endif
                d.Vc(PRS,k,j,i) = d.Vc(RHO,k,j,i)/cs2*gammaIdeal;
                d.Vc(VX1,k,j,i) = 1e-4*(idfx::randm()-0.5);
                d.Vc(VX2,k,j,i) = shear*x;
                d.Vc(VX3,k,j,i) = 1e-4*(idfx::randm()-0.5);

                d.Vs(BX1s,k,j,i) = 0;
                d.Vs(BX2s,k,j,i) = B0y;
                d.Vs(BX3s,k,j,i) = B0z;

            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}


// Analyse data to produce an output

void MakeAnalysis(DataBlock & data) {

}



// Do a specifically designed user step in the middle of the integration
void ComputeUserStep(DataBlock &data, real t, real dt) {

}
