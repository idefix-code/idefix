#include "idefix.hpp"
#include "setup.hpp"
#include "column.hpp"
#include "lookupTable.hpp"

real R0Glob;
real h0Glob;
real hpowGlob;
real rho0Glob;
real rhominGlob;
real T0Glob;
real muGlob;
real gammaGlob;
real rsGlob;
real TsGlob;
real kappaGlob;
real kappairrGlob;
std::string kappairrtypeGlob;
std::string kappatypeGlob;
std::string xitypeGlob;

Column *columnGlob;
Column *columnGlob2;
LookupTable<1> *kappairrtableGlob;

void MyKappa(DataBlock &data, IdefixArray3D<real> &kappap, IdefixArray3D<real> &kappar) {

  idefix_for("MyKappa",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                kappap(k,j,i) = 1.e1;
                kappar(k,j,i) = 1.e1;
              });
}

void MyXi(DataBlock &data, IdefixArray3D<real> &xi) {

  idefix_for("MyXi",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                xi(k,j,i) = 0.;
              });
}


void FluxBoundary(Fluid<DefaultPhysics> *hydro, int dir, BoundarySide side, const real t) {
    idfx::pushRegion("FluxInternal");

    auto *data = hydro->data;

    IdefixArray4D<real> Flux = data->hydro->FluxRiemann[dir];

    if( dir== IDIR ) {
      idefix_for("FluxInternal",
                  0, data->np_tot[KDIR],
                  0, data->np_tot[JDIR],
                  0, data->np_tot[IDIR],
         KOKKOS_LAMBDA (int k, int j, int i) {
           Flux(RHO, k, j, i) = 0.0;
           Flux(MX1, k, j, i) = 0.0;
           Flux(MX2, k, j, i) = 0.0;
           Flux(MX3, k, j, i) = 0.0;
           Flux(ENG, k, j, i) = 0.0;
      });
    }

    if( dir== JDIR ) {
      idefix_for("FluxInternal",
                  0, data->np_tot[KDIR],
                  0, data->np_tot[JDIR],
                  0, data->np_tot[IDIR],
         KOKKOS_LAMBDA (int k, int j, int i) {
           Flux(RHO, k, j, i) = 0.0;
           Flux(MX1, k, j, i) = 0.0;
           Flux(MX2, k, j, i) = 0.0;
           Flux(MX3, k, j, i) = 0.0;
           Flux(ENG, k, j, i) = 0.0;
      });
    }

    if( dir== KDIR ) {
      idefix_for("FluxInternal",
                  0, data->np_tot[KDIR],
                  0, data->np_tot[JDIR],
                  0, data->np_tot[IDIR],
         KOKKOS_LAMBDA (int k, int j, int i) {
           Flux(RHO, k, j, i) = 0.0;
           Flux(MX1, k, j, i) = 0.0;
           Flux(MX2, k, j, i) = 0.0;
           Flux(MX3, k, j, i) = 0.0;
           Flux(ENG, k, j, i) = 0.0;
      });
    }

    idfx::popRegion();

}

void InternalBoundary(Hydro *hydro, const real t) {
  auto *data = hydro->data;
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray4D<real> Uc = hydro->Uc;
  IdefixArray1D<real> x1=data->x[IDIR];
  IdefixArray1D<real> x2=data->x[JDIR];

  real R0 = R0Glob;
  real h0 = h0Glob;
  real hpow = hpowGlob;
  real rho0 = rho0Glob;
  real rhomin = rhominGlob;
  real unit_density = idfx::units.GetDensity();


    idefix_for("InternalBoundary",
      0, data->np_tot[KDIR],
      0, data->np_tot[JDIR],
      0, data->np_tot[IDIR],
            KOKKOS_LAMBDA (int k, int j, int i) {

              real R = FMAX(x1(i)*std::sin(x2(j)),1.);
              real z2 = std::pow(x1(i)*std::cos(x2(j)),2.);
              real H = h0*std::pow(R/R0,hpow);

              Vc(RHO,k,j,i) = (rho0*(R0/R)*std::exp(-0.25*M_PI*z2/(H*H))+rhomin)/unit_density;
              Vc(VX1,k,j,i) = 0.;
              Vc(VX2,k,j,i) = 0.;
              Vc(VX3,k,j,i) = 0.;

              Uc(RHO,k,j,i) = (rho0*(R0/R)*std::exp(-0.25*M_PI*z2/(H*H))+rhomin)/unit_density;
              Uc(VX1,k,j,i) = 0.;
              Uc(VX2,k,j,i) = 0.;
              Uc(VX3,k,j,i) = 0.;
            });

}

// Compute user variables which will be written in vtk files
void ComputeUserVars(DataBlock & data, UserDefVariablesContainer &variables) {
  // Mirror data on Host
  DataBlockHost d(data);

  // Sync it
  d.SyncFromDevice();

  IdefixHostArray1D<real> x1=d.x[IDIR];
  IdefixHostArray1D<real> x1l=d.xl[IDIR];
  IdefixHostArray1D<real> x2=d.x[JDIR];
  IdefixHostArray1D<real> x2l=d.xl[JDIR];
  IdefixHostArray3D<real> A1=d.A[IDIR];
  IdefixHostArray3D<real> A2=d.A[JDIR];
  IdefixHostArray3D<real> dV=d.dV;

  IdefixArray3D<real> tau;
  IdefixArray4D<real> Vc=(&data)->hydro->Vc;

  IdefixHostArray3D<real> divF  = variables["divF"];
  IdefixHostArray3D<real> A1_out  = variables["A1"];

  std::string kappairrType = kappairrtypeGlob;

  real rs = rsGlob;
  real Ts = TsGlob;
  real flux_pre = std::pow(rs/idfx::units.GetLength(),2.)*idfx::units.sigma_sb*std::pow(Ts,4.)/idfx::units.GetLength();

  // Make references to the user-defined arrays (variables is a container of IdefixHostArray3D)
  // Note that the labels should match the variable names in the input file

  columnGlob->ComputeColumn(Vc,RHO);
  tau = columnGlob->GetColumn();

  IdefixHostArray3D<real> var_tau = variables["tau"];
  IdefixHostArray3D<real> var_dV = variables["dV"];

  Kokkos::deep_copy(var_tau, tau);
  Kokkos::deep_copy(var_dV, dV);

  IdefixArray3D<real> rho("rho",d.np_tot[KDIR],d.np_tot[JDIR],d.np_tot[IDIR]);
  idefix_for("init rho",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
    KOKKOS_LAMBDA(int k, int j, int i) {
      rho(k,j,i) = Vc(RHO,k,j,i);
    });

  columnGlob2->ComputeColumn(rho);

  if(kappairrType=="constant") {

    real kappa_irr = kappairrGlob;
    real kirr = kappa_irr*idfx::units.GetDensity()*idfx::units.GetLength();

    for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
          A1_out(k,j,i) = A1(k,j,i);

          real Fim = std::exp(-kirr*var_tau(k,j,i-1))/std::pow(x1l(i),2.);
          real Fip = std::exp(-kirr*var_tau(k,j,i))/std::pow(x1l(i+1),2.);

          divF(k,j,i) = (Fip*A1(k,j,i+1)-Fim*A1(k,j,i));
          divF(k,j,i) *= flux_pre;
          divF(k,j,i) /= dV(k,j,i);

        }
      }
    }
  } else if (kappairrType=="usertable") {

    for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
          A1_out(k,j,i) = A1(k,j,i);

          real logtaum = std::log10(FMAX(var_tau(k,j,i-1)*idfx::units.GetDensity()*idfx::units.GetLength(),1.e-15));
          real Fim = pow(10.,kappairrtableGlob->GetHost(&logtaum))*A1(k,j,i)/std::pow(x1l(i),2.);
          real logtaup = std::log10(FMAX(var_tau(k,j,i)*idfx::units.GetDensity()*idfx::units.GetLength(),1.e-15));
          real Fip = pow(10.,kappairrtableGlob->GetHost(&logtaup))*A1(k,j,i+1)/std::pow(x1l(i+1),2.);
          divF(k,j,i)  = flux_pre*(Fip-Fim)/dV(k,j,i);

        }
      }
    }
  }
  Kokkos::deep_copy(variables["divF"], divF);
  Kokkos::deep_copy(variables["A1"], A1_out);

}


// Default constructor
// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output)
{
  DataBlockHost d(data);
  auto Vc = data.hydro->Vc;

  // Mirror data on Host
  R0Glob = input.Get<real>("Setup","R0",0);
  h0Glob = input.Get<real>("Setup","h0",0);
  hpowGlob = input.Get<real>("Setup","hpow",0);
  rho0Glob = input.Get<real>("Setup","rho0",0);
  rhominGlob = input.Get<real>("Setup","rhomin",0);
  T0Glob = input.Get<real>("Setup","T0",0);
  muGlob = input.Get<real>("Hydro","mu",0);
  gammaGlob=data.hydro->eos->GetGamma();
  rsGlob=input.Get<real>("Rad","irr",1);
  TsGlob=input.Get<real>("Rad","irr",2);

  kappairrtypeGlob = input.Get<std::string>("Rad","irr",0);
  kappatypeGlob = input.Get<std::string>("Rad","kappa",0);
  xitypeGlob = input.Get<std::string>("Rad","xi",0);

  if (kappairrtypeGlob == "constant") {
    kappairrGlob = input.Get<real>("Rad","irr",3);
  } else if (kappairrtypeGlob == "usertable") {
    std::string irr_file = input.Get<std::string>("Rad","irr",4);
    kappairrtableGlob = new LookupTable<1>(irr_file,',');
  }

  columnGlob = new Column(IDIR,1,&data);
  columnGlob2 = new Column(IDIR,1,&data);

  if (kappatypeGlob == "userfunc") {
    data.radiation[0]->EnrollKappa(&MyKappa);
  }
  if (xitypeGlob == "userfunc") {
    data.radiation[0]->EnrollXi(&MyXi);
  }

  data.hydro->EnrollInternalBoundary(&InternalBoundary);
  data.hydro->EnrollFluxBoundary(&FluxBoundary);
  output.EnrollUserDefVariables(&ComputeUserVars);
}

Setup::~Setup() {
  delete columnGlob;
  delete columnGlob2;
}
// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);

    real R0 = R0Glob;
    real h0 = h0Glob;
    real hpow = hpowGlob;
    real rho0 = rho0Glob;
    real rhomin = rhominGlob;
    real T0 = T0Glob;
    real mu = muGlob;


    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {

              real R = FMAX(d.x[IDIR](i)*std::sin(d.x[JDIR](j)),1.);
              real z2 = std::pow(d.x[IDIR](i)*std::cos(d.x[JDIR](j)),2.);
              real H = h0*std::pow(R/R0,hpow);
              real T = T0;

              d.Vc(RHO,k,j,i) = (rho0*(R0/R)*std::exp(-0.25*M_PI*z2/(H*H))+rhomin)/idfx::units.GetDensity();
              d.Vc(PRS,k,j,i) = d.Vc(RHO,k,j,i)*idfx::units.GetDensity()*idfx::units.k_B*T/(mu*idfx::units.u)/idfx::units.GetEnergy();
              d.Vc(VX1,k,j,i) = 0.;
              d.Vc(VX2,k,j,i) = 0.;
              d.Vc(VX3,k,j,i) = 0.;
              d.RadVc[0](ER,k,j,i) = idfx::units.ar*std::pow(T,4.)/idfx::units.GetEnergy();
              d.RadVc[0](FR1,k,j,i) = ZERO_F;
              d.RadVc[0](FR2,k,j,i) = ZERO_F;
              d.RadVc[0](FR3,k,j,i) = ZERO_F;

              d.RadUc[0](ER,k,j,i) = idfx::units.ar*std::pow(T,4.)/idfx::units.GetEnergy();
              d.RadUc[0](FR1,k,j,i) = ZERO_F;
              d.RadUc[0](FR2,k,j,i) = ZERO_F;
              d.RadUc[0](FR3,k,j,i) = ZERO_F;
            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}
