#include "idefix.hpp"
#include "setup.hpp"

real epsilonGlob;
real betaGlob;
real HidealGlob;
real AmMidGlob;
real gammaGlob;
real densityFloorGlob;
real alphaGlob;


void MySoundSpeed(DataBlock &data, const real t, IdefixArray3D<real> &cs) {
  IdefixArray1D<real> r=data.x[IDIR];
  IdefixArray1D<real> th=data.x[JDIR];
  real epsilon = epsilonGlob;
  idefix_for("MySoundSpeed",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real R = r(i)*sin(th(j));
                cs(k,j,i) = epsilon/sqrt(R);
              });
}

void MyViscosity(DataBlock &data, const real t, IdefixArray3D<real> &eta1, IdefixArray3D<real> &eta2) {
  IdefixArray4D<real> Vc=data.hydro->Vc;
  IdefixArray1D<real> r=data.x[IDIR];
  IdefixArray1D<real> th=data.x[JDIR];
  real epsilon = epsilonGlob;
  real alpha = alphaGlob;
  idefix_for("MyViscosity",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real R = r(i)*sin(th(j));
                real cs = epsilon/sqrt(R);
                eta1(k,j,i) = alpha*cs*epsilon*R*Vc(RHO,k,j,i);
                eta2(k,j,i) = ZERO_F;
              });

}

void FargoVelocity(DataBlock &data, IdefixArray2D<real> &Vphi) {
  IdefixArray1D<real> x1 = data.x[IDIR];

  idefix_for("FargoVphi",0,data.np_tot[KDIR], 0, data.np_tot[IDIR],
      KOKKOS_LAMBDA (int k, int i) {
      Vphi(k,i) = 1.0/sqrt(x1(i));
  });
}

// User-defined boundaries
void UserdefBoundary(Hydro *hydro, int dir, BoundarySide side, real t) {
  auto *data = hydro->data;
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  real epsilon=epsilonGlob;
    if(dir==IDIR) {


        int ighost;
        if(side == left) {
            ighost = data->beg[IDIR];
        }
        else if(side==right) {
            ighost = data->end[IDIR]-1;
        }

        hydro->boundary->BoundaryFor("UserDefBoundary", dir, side,
          KOKKOS_LAMBDA (int k, int j, int i) {
              real R=x1(i)*sin(x2(j));
              real z=x1(i)*cos(x2(j));
              real Vk = 1.0/sqrt(R);
              real cs2=(epsilon*Vk)*(epsilon*Vk);

              Vc(RHO,k,j,i) = 1.0/(R*sqrt(R))*exp(1.0/(cs2)*(1/x1(i)-1/R));
              Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost);
              Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost);
              Vc(VX3,k,j,i) = Vk*sqrt(R/x1(i)-2.5*epsilon*epsilon);
          });
    }

    if( dir==JDIR) {
        int jghost;
        if(side == left) {
            jghost = data->beg[JDIR];
        }
        else if(side==right) {
            jghost = data->end[JDIR]-1;
        }


        hydro->boundary->BoundaryFor("UserDefBoundary", dir, side,
            KOKKOS_LAMBDA (int k, int j, int i) {
              real R=x1(i)*sin(x2(j));
              real z=x1(i)*cos(x2(j));
              real Vk = 1.0/sqrt(R);
              real cs2=(epsilon*Vk)*(epsilon*Vk);

              Vc(RHO,k,j,i) = 1.0/(R*sqrt(R))*exp(1.0/(cs2)*(1/x1(i)-1/R));
              Vc(VX1,k,j,i) = Vc(VX1,k,jghost,i);
              Vc(VX2,k,j,i) = ZERO_F;
              Vc(VX3,k,j,i) = Vk*sqrt(R/x1(i)-2.5*epsilon*epsilon);

            });
    }

}

// Default constructor


// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  // Set the function for userdefboundary
  data.hydro->EnrollUserDefBoundary(&UserdefBoundary);
  data.hydro->EnrollIsoSoundSpeed(&MySoundSpeed);
  data.hydro->viscosity->EnrollViscousDiffusivity(&MyViscosity);
  if(data.haveFargo)
    data.fargo->EnrollVelocity(&FargoVelocity);
  epsilonGlob = input.Get<real>("Setup","epsilon",0);
  alphaGlob = input.Get<real>("Setup","alpha",0);
  idfx::cout << "alpha= " << alphaGlob << std::endl;
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);
    real r,th;
    real epsilon=epsilonGlob;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                r=d.x[IDIR](i);
                th=d.x[JDIR](j);
                real R=r*sin(th);
                real z=r*cos(th);
                real Vk=1.0/sqrt(R);

                real cs2=(epsilon*Vk)*(epsilon*Vk);

                d.Vc(RHO,k,j,i) = 1.0/(R*sqrt(R))*exp(1.0/(cs2)*(1/r-1/R));
                d.Vc(VX1,k,j,i) = 0.0;
                d.Vc(VX2,k,j,i) = 0.0;
                d.Vc(VX3,k,j,i) = Vk*sqrt(R/r-2.5*epsilon*epsilon);
            }
        }
    }


    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}
