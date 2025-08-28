#include "idefix.hpp"
#include "setup.hpp"

real udenGlob;
real ulenGlob;
real uvelGlob;
real gammaGlob;
real cs_vescGlob;
real va_vescGlob;
real k0_Glob;
real ParkerWind(real);

const real kB{1.3807e-16};
const real mp{1.6726e-24};

void MyClessThermalConductivity(DataBlock &data, const real t, std::vector<IdefixArray3D<real>> &userdefArr) {
  IdefixArray4D<real> Vc = data.hydro->Vc;
  IdefixArray1D<real> x1 = data.x[IDIR];
  IdefixArray1D<real> x2 = data.x[JDIR];

  IdefixArray3D<real> kparArr = userdefArr.at(0);
  IdefixArray3D<real> knorArr = userdefArr.at(1);
  IdefixArray3D<real> clessAlpha = userdefArr.at(2);
  IdefixArray3D<real> clessBeta  = userdefArr.at(3);

  real norm = mp*0.5/(udenGlob*uvelGlob*ulenGlob*kB);
  real uTemp=0.5*uvelGlob*uvelGlob*mp/kB;
  real k0 = k0_Glob*norm;
  idefix_for("MyClessThConductivity",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
         KOKKOS_LAMBDA (int k, int j, int i) {
           kparArr(k,j,i) = k0*pow(Vc(PRS,k,j,i)/Vc(RHO,k,j,i)*uTemp,2.5);
           knorArr(k,j,i) = 0.;
           clessAlpha(k,j,i) = (1.0-tanh(x1(i)-10))/2;
           clessBeta(k,j,i)  = -1.5;
         });
}

// User-defined boundaries
void UserDefBoundary(Hydro *hydro, int dir, BoundarySide side, real t) {

  DataBlock *data = hydro->data;

  if( (dir==IDIR) && (side == left)) {
    IdefixArray4D<real> Vc = hydro->Vc;
    IdefixArray4D<real> Vs = hydro->Vs;
    IdefixArray1D<real> x1 = data->x[IDIR];

    real rc,vwind0;
    real cs=cs_vescGlob*sqrt(2.);
    real va_vesc = va_vescGlob;
    real mu = va_vesc * sqrt(2.);
    real PonRho;

    PonRho = cs*cs;
    rc = 0.25 / (cs_vescGlob*cs_vescGlob);
    vwind0 = ParkerWind(1./rc) * cs;

    hydro->boundary->BoundaryFor("UserDefBoundary", dir, side,
                 KOKKOS_LAMBDA (int k, int j, int i) {
                   real r = x1(i);

                   Vc(RHO,k,j,i) = vwind0/(vwind0 * r * r);
                   Vc(PRS,k,j,i) = PonRho * Vc(RHO, k, j, i);
                   Vc(VX1,k,j,i) = vwind0;
                   Vc(VX2,k,j,i) = 0.0;
                   Vc(VX3,k,j,i) = 0.0;
                   Vc(BX1,k,j,i) = mu / (r*r);
                   Vc(BX2,k,j,i) = 0.0;
                   Vc(BX3,k,j,i) = 0.0;

                 });
  }
}

// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  // Set the function for userdefboundary
  data.hydro->EnrollUserDefBoundary(&UserDefBoundary);
  data.hydro->bragThermalDiffusion->EnrollBragThermalDiffusivity(&MyClessThermalConductivity);
  gammaGlob=input.Get<real>("Hydro", "gamma", 0);
  udenGlob=input.Get<real>("Setup", "UNIT_DENSITY",0);
  ulenGlob=input.Get<real>("Setup", "UNIT_LENGTH",0);
  uvelGlob=input.Get<real>("Setup", "UNIT_VELOCITY",0);
  cs_vescGlob=input.Get<real>("Setup", "cs_vesc", 0);
  va_vescGlob=input.Get<real>("Setup", "va_vesc", 0);
  k0_Glob = input.Get<real>("Setup", "k0", 0);
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
  // Create a host copy
  DataBlockHost d(data);

  real r,rl;
  real PonRho, vwind0, rc;
  real cs=cs_vescGlob*sqrt(2.);


  rc = 0.25 / (cs_vescGlob*cs_vescGlob);
  vwind0 = ParkerWind(1./rc) * cs;
  PonRho = cs*cs;
  real mu = va_vescGlob * sqrt(2.);

  for(int k = 0; k < d.np_tot[KDIR] ; k++) {
    for(int j = 0; j < d.np_tot[JDIR] ; j++) {
      for(int i = 0; i < d.np_tot[IDIR] ; i++) {
        r=d.x[IDIR](i);

        real vwind;

        vwind = ParkerWind(r/rc) * cs;

        d.Vc(RHO,k,j,i) = 1.0*vwind0/(vwind * r * r);
        d.Vc(PRS,k,j,i) = PonRho * d.Vc(RHO, k, j, i);
        d.Vc(VX1,k,j,i) = vwind;
        d.Vc(VX2,k,j,i) = 0.0;
        d.Vc(VX3,k,j,i) = 0.0;

        rl=d.xl[IDIR](i); // Radius on the left side of the cell
        d.Vs(BX1s, k, j, i) = mu / (rl*rl);

      }
    }
  }
  // Send it all, if needed
  d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {
}

/**************************************************/
real ParkerWind(real x)
/*  Parker wind velocity in unit of iso sound speed
    x = radius / critical radius.
**************************************************/
{
  real v, f;
  real vref;

  v = 1e-7;
  f = v*v-2*log(v)-4/x-4*log(x)+3;
  if (x>1) {v=10;}
  while (fabs(f) > 1e-10){
    vref = v;
    v = v - 0.5*f/(v-1/v);
    while (v < 0){
      v = (v + vref)/2;
    }
    f = v*v-2*log(v)-4/x-4*log(x)+3;
  }
  return v;
}
