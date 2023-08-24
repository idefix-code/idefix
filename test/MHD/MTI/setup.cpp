#include "idefix.hpp"
#include "setup.hpp"
#include "analysis.hpp"

// Definition of the constants and parameters of the problem
static real T0 = 1.;
static real H = 3.;
static real L = 0.1;
static real g0 = 1.;

static real beta0 = 1e12;
static real rho0 = 1.; //initial density
static real P0 = rho0*T0; //initial pressure
static real B0 = 100.*sqrt(1./beta0);

static real vth0 = 1.;

static real n = 1.;
static real kn = 2.*n*M_PI/L;
static real ksiGlob;
static real prGlob;

Analysis *analysis;

void AnalysisFunction(DataBlock &data) {
  analysis->PerformAnalysis(data);
}

// Implement uniform gravitationnal field
void Potential(DataBlock& data, const real t, IdefixArray1D<real>& x1, IdefixArray1D<real>& x2, IdefixArray1D<real>& x3, IdefixArray3D<real>& phi) {
  real g0 = 1.;
  idefix_for("Potential", 0, data.np_tot[KDIR], 0, data.np_tot[JDIR], 0, data.np_tot[IDIR],
      KOKKOS_LAMBDA (int k, int j, int i) {
          phi(k,j,i) = g0*x2(j);
  });
}


void MyBragThermalConductivity(DataBlock &data, const real t, IdefixArray3D<real> &kparArr, IdefixArray3D<real> &knorArr) {
  IdefixArray4D<real> Vc = data.hydro->Vc;
  IdefixArray1D<real> x2 = data.x[JDIR];
  real ksi = ksiGlob;
  idefix_for("MyThConductivity",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                kparArr(k,j,i) = ksi*Vc(RHO,k,j,i);
                knorArr(k,j,i) = 0.;
  });
}

void MyViscosity(DataBlock &data, const real t, IdefixArray3D<real> &etaBrag) {
  IdefixArray1D<real> x1 = data.x[IDIR];
  IdefixArray1D<real> x2 = data.x[JDIR];
  real ksi = ksiGlob;
  real pr = prGlob;
  IdefixArray4D<real> Vc = data.hydro->Vc;
  idefix_for("MyViscosity",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                etaBrag(k,j,i) = ksi*pr*Vc(RHO,k,j,i);
              });
}

// Define our own boundary conditions. Basically, homogene Neummann (i.e. symmetry) on polar and azimuthal velocity, Dirichlet on the thermo fields.
void UserDefBoundary(Hydro *hydro, int dir, BoundarySide side, real t) {
  auto *data = hydro->data;
  real T0 = 1.;
  real H = 3.;
  real g0 = 1.;
  real rho0 = 1.; //initial density
  real P0 = rho0*T0; //initial pressure
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray4D<real> Vs = hydro->Vs;
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  int jbeg = (side == left) ? 0 : data->end[JDIR];
  int jend = (side == left) ? data->beg[JDIR] : data->np_tot[JDIR];
  hydro->boundary->BoundaryFor("UserDefBoundary", dir, side,
      KOKKOS_LAMBDA (int k, int j, int i) {
          int jref = (side == left) ? 2*jend - j -1 : jbeg + jend -j - 3;
          Vc(PRS,k,j,i) = P0*pow(1.-x2(j)/H, H*g0/T0);
          Vc(RHO,k,j,i) =rho0*pow(1.-x2(j)/H, H*g0/T0-1.) - (Vc(RHO,k,jref,i) - rho0*pow(1.-x2(jref)/H, H*g0/T0-1.));

          Vc(VX1,k,j,i) = Vc(VX1,k,jref,i);
          Vc(VX2,k,j,i) = -Vc(VX2,k,jref,i);
      });

  hydro->boundary->BoundaryForX2s("UserDefBoundaryBX2s", dir, side,
      KOKKOS_LAMBDA (int k, int j, int i) {
          int jref = (side == left) ? 2*jend - j -1 : jbeg + jend -j - 3;
          Vs(BX2s,k,j,i) = Vs(BX2s,k,jref,i);
  });
}

// Initialisation routine. Can be used to allocate Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  analysis = new Analysis(input, grid, data, output,std::string("average.dat"));
  output.EnrollAnalysis(&AnalysisFunction);
  // Reset analysis if required
  if(!input.restartRequested) {
    analysis->ResetAnalysis();
  }

  data.gravity->EnrollPotential(&Potential);
  data.hydro->EnrollUserDefBoundary(&UserDefBoundary);
  data.hydro->bragThermalDiffusion->EnrollBragThermalDiffusivity(&MyBragThermalConductivity);
  data.hydro->bragViscosity->EnrollBragViscousDiffusivity(&MyViscosity);
  ksiGlob = input.Get<real>("Setup","ksi",0);
  prGlob = input.Get<real>("Setup","pr",0);
}

// This routine initializes the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
  // Create a host copy
  DataBlockHost d(data);

  for(int k = 0; k < d.np_tot[KDIR] ; k++) {
    for(int j = 0; j < d.np_tot[JDIR] ; j++) {
      real x2 = d.x[JDIR](j);
      real x2s = d.xl[JDIR](j);
      for(int i = 0; i < d.np_tot[IDIR] ; i++) {
        real x1=d.x[IDIR](i);
        real x1s=d.xl[IDIR](i);

        d.Vc(PRS,k,j,i) = P0*pow(1.-x2/H, H*g0/T0);
        d.Vc(RHO,k,j,i) = rho0*pow(1.-x2/H, H*g0/T0-1.);

        d.Vc(VX1,k,j,i) = -1e-4*SIN(kn*x1)*COS(kn*x2)*vth0;
        d.Vc(VX2,k,j,i) = 1e-4*COS(kn*x1)*SIN(kn*x2)*vth0;

        d.Vs(BX1s,k,j,i) = B0;
        d.Vs(BX2s,k,j,i) = 0.;
      }
    }
  }
  d.SyncToDevice();
}
