#include "idefix.hpp"
#include "setup.hpp"
#include "analysis.hpp"
#include "dumpImage.hpp"

// Definition of the constants and parameters of the problem
const real T0 = 1.;
const real H = 3.;
const real L = 0.1;
const real g0 = 1.;

const real beta0 = 1e12;
const real rho0 = 1.; //initial density
const real P0 = rho0*T0; //initial pressure
const real B0 = 100.*sqrt(1./beta0);

const real vth0 = 1.;

const real n = 1.;
const real kn = 2.*n*M_PI/L;
static real ksiGlob;
static real prGlob;

Analysis *analysis;

bool fromDump;

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

  hydro->boundary->BoundaryForX1s("UserDefBoundaryBX1s", dir, side,
      KOKKOS_LAMBDA (int k, int j, int i) {
          //int jref = (side == left) ? 2*jend - j -1 : jbeg + jend -j - 3;
          int jref = (side == left) ? 2*jend - j -1 : 2*jbeg -j -1;
          Vs(BX1s,k,j,i) = Vs(BX1s,k,jref,i);
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
  fromDump = input.GetOrSet<bool>("Setup","fromDump",0,false);  // Whether we build our own initial condition or we use the pre-defined one
}

// This routine initializes the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
  // Create a host copy
  DataBlockHost d(data);

  if(fromDump) {
    DumpImage image("init.dmp",&data);

    for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {

          // Note that the restart dump array only contains the full (global) active domain
          // (i.e. it excludes the boundaries, but it is not decomposed accross MPI procs)
          int iglob=i-2*d.beg[IDIR]+d.gbeg[IDIR];
          int jglob=j-2*d.beg[JDIR]+d.gbeg[JDIR];
          int kglob=k-2*d.beg[KDIR]+d.gbeg[KDIR];

          d.Vc(RHO,k,j,i) = image.arrays["Vc-RHO"](kglob,jglob,iglob);
          d.Vc(PRS,k,j,i) = image.arrays["Vc-PRS"](kglob,jglob,iglob);
          d.Vc(VX1,k,j,i) = image.arrays["Vc-VX1"](kglob,jglob,iglob);
          d.Vc(VX2,k,j,i) = image.arrays["Vc-VX2"](kglob,jglob,iglob);
        }}}
    for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR]+1 ; i++) {
          int iglob=i-2*d.beg[IDIR]+d.gbeg[IDIR];
          int jglob=j-2*d.beg[JDIR]+d.gbeg[JDIR];
          int kglob=k-2*d.beg[KDIR]+d.gbeg[KDIR];
          d.Vs(BX1s,k,j,i) = image.arrays["Vs-BX1s"](kglob,jglob,iglob);
    }}}
    for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR]+1 ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
          int iglob=i-2*d.beg[IDIR]+d.gbeg[IDIR];
          int jglob=j-2*d.beg[JDIR]+d.gbeg[JDIR];
          int kglob=k-2*d.beg[KDIR]+d.gbeg[KDIR];
          d.Vs(BX2s,k,j,i) = image.arrays["Vs-BX2s"](kglob,jglob,iglob);
    }}}
  } else {
    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
      for(int j = 0; j < d.np_tot[JDIR] ; j++) {
        real x2 = d.x[JDIR](j);
        real x2s = d.xl[JDIR](j);
        for(int i = 0; i < d.np_tot[IDIR] ; i++) {
          real x1=d.x[IDIR](i);
          real x1s=d.xl[IDIR](i);

          d.Vc(PRS,k,j,i) = P0*pow(1.-x2/H, H*g0/T0);
          d.Vc(RHO,k,j,i) = rho0*pow(1.-x2/H, H*g0/T0-1.);

          d.Vc(VX1,k,j,i) = -1e-4*std::sin(kn*x1)*std::cos(kn*x2)*vth0;
          d.Vc(VX2,k,j,i) = 1e-4*std::cos(kn*x1)*std::sin(kn*x2)*vth0;

          d.Vs(BX1s,k,j,i) = B0;
          d.Vs(BX2s,k,j,i) = 0.;
        }
      }
    }
  }
  d.SyncToDevice();
}
