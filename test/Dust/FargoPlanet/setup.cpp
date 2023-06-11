#include <algorithm>
#include "idefix.hpp"
#include "setup.hpp"

real sigmaSlopeGlob;
real h0Glob;
real HidealGlob;
real gammaGlob;
real densityFloorGlob;
real alphaGlob;


void MySoundSpeed(DataBlock &data, const real t, IdefixArray3D<real> &cs) {
  IdefixArray1D<real> x1=data.x[IDIR];
  real h0 = h0Glob;
  idefix_for("MySoundSpeed",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real R = x1(i);
                cs(k,j,i) = h0/sqrt(R);
              });
}

// drag coefficient assuming vertical hydrostatic equilirbium
// And using the fact that Vc(RHO) is the surface density
// This assumes that tau_s = beta * sqrt(2pi) h / (cs*Sigma)
// So that the Stokes number St = Omega tau_s = sqrt(2pi) * beta / Sigma

// In this setup, we assume rho_mid(R=1)=1, so that Sigma(R=1)=sqrt(2pi)*h0
// hence beta = h0 * St(R=1)

// Assuming Epstein drag, the particle size is related to St through:
// a = 20 cm *  St * (Sigma/100 g.cm^-2) / (rho_s / 2 g.cm^-3)

void MyDrag(DataBlock *data, real beta, IdefixArray3D<real> &gamma) {
  // Compute the drag coefficient gamma from the input beta
  real h0 = h0Glob;
  auto x1 = data->x[IDIR];

  idefix_for("MyDrag",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      real R = x1(i);
      real h = h0*R;
      real cs = h0/sqrt(R);


      gamma(k,j,i) = cs/(beta*sqrt(2.0*M_PI)*h);

    });
}

// User-defined boundaries
void UserdefBoundary(Hydro *hydro, int dir, BoundarySide side, real t) {
  IdefixArray4D<real> Vc = hydro->Vc;
  auto *data = hydro->data;
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x3 = data->x[KDIR];
  if(dir==IDIR) {
    int ighost,ibeg,iend;
    if(side == left) {
      ighost = data->beg[IDIR];
      ibeg = 0;
      iend = data->beg[IDIR];
      idefix_for("UserDefBoundary",0,data->np_tot[KDIR],0,data->np_tot[JDIR],ibeg,iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          real R=x1(i);
          real z=x3(k);
          real Vk = 1.0/sqrt(R);

          Vc(RHO,k,j,i) = Vc(RHO,k,j,2*ighost - i +1);
          Vc(VX1,k,j,i) = - Vc(VX1,k,j,2*ighost - i +1);
          Vc(VX2,k,j,i) = Vk;
          Vc(VX3,k,j,i) = Vc(VX3,k,j,2*ighost - i +1);
        });
    }
    else if(side==right) {
      ighost = data->end[IDIR]-1;
      ibeg=data->end[IDIR];
      iend=data->np_tot[IDIR];
      idefix_for("UserDefBoundary",0,data->np_tot[KDIR],0,data->np_tot[JDIR],ibeg,iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          real R=x1(i);
          real z=x3(k);
          real Vk = 1.0/sqrt(R);

          Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost);
          Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost);
          Vc(VX2,k,j,i) = Vk;
          Vc(VX3,k,j,i) = Vc(VX3,k,j,ighost);
        });
    }
  }
}

void UserdefBoundaryDust(Fluid<DustPhysics> *dust, int dir, BoundarySide side, real t) {
  IdefixArray4D<real> Vc = dust->Vc;
  auto data = dust->data;
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x3 = data->x[KDIR];
  if(dir==IDIR) {
    int ighost,ibeg,iend;
    if(side == left) {
      ighost = data->beg[IDIR];
      ibeg = 0;
      iend = data->beg[IDIR];
      idefix_for("UserDefBoundary",0,data->np_tot[KDIR],0,data->np_tot[JDIR],ibeg,iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          real R=x1(i);
          real z=x3(k);
          real Vk = 1.0/sqrt(R);

          Vc(RHO,k,j,i) = Vc(RHO,k,j,2*ighost - i +1);
          Vc(VX1,k,j,i) = - Vc(VX1,k,j,2*ighost - i +1);
          Vc(VX2,k,j,i) = Vk;
          Vc(VX3,k,j,i) = Vc(VX3,k,j,2*ighost - i +1);
        });
    }
    else if(side==right) {
      ighost = data->end[IDIR]-1;
      ibeg=data->end[IDIR];
      iend=data->np_tot[IDIR];
      idefix_for("UserDefBoundary",0,data->np_tot[KDIR],0,data->np_tot[JDIR],ibeg,iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          real R=x1(i);
          real z=x3(k);
          real Vk = 1.0/sqrt(R);

          Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost);
          Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost);
          Vc(VX2,k,j,i) = Vk;
          Vc(VX3,k,j,i) = Vc(VX3,k,j,ighost);
        });
    }
  }
}



void FargoVelocity(DataBlock &data, IdefixArray2D<real> &Vphi) {
  IdefixArray1D<real> x1 = data.x[IDIR];

  idefix_for("FargoVphi",0,data.np_tot[KDIR], 0, data.np_tot[IDIR],
      KOKKOS_LAMBDA (int k, int i) {
      Vphi(k,i) = 1.0/sqrt(x1(i));
  });
}


// Default constructor
// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output)// : m_planet(0)//, Planet &planet)
{
  // Set the function for userdefboundary
  data.hydro->EnrollUserDefBoundary(&UserdefBoundary);
  if(data.haveDust) {
    int nSpecies = data.dust.size();
    for(int n = 0 ; n < nSpecies ; n++) {
      data.dust[n]->EnrollUserDefBoundary(&UserdefBoundaryDust);
      data.dust[n]->drag->EnrollUserDrag(&MyDrag);
    }
  }

  data.hydro->EnrollIsoSoundSpeed(&MySoundSpeed);
  if(data.haveFargo)
    data.fargo->EnrollVelocity(&FargoVelocity);
  sigmaSlopeGlob = input.Get<real>("Setup","sigmaSlope",0);
  h0Glob = input.Get<real>("Setup","h0",0);
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);
    real h0=h0Glob;
    real sigmaSlope = sigmaSlopeGlob;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                real R=d.x[IDIR](i);
                real z=d.x[KDIR](k);
                real Vk=1.0/sqrt(R);

                real cs2=(h0*Vk)*(h0*Vk);

                d.Vc(RHO,k,j,i) = pow(R,-sigmaSlope-1) * exp(1.0/ (cs2) * (1.0/sqrt(R*R+z*z)-1.0/R)) ;
                d.Vc(VX1,k,j,i) = 0.0;//+sin(8*d.x[JDIR](j));
                d.Vc(VX3,k,j,i) = 0.0;
                d.Vc(VX2,k,j,i) = Vk*sqrt( R / sqrt(R*R + z*z) -(2.0+sigmaSlope)*h0*h0);

                for(int n = 0 ; n < data.dust.size() ; n++) {
                  d.dustVc[n](RHO,k,j,i) = 1e-2*d.Vc(RHO,k,j,i);
                  d.dustVc[n](VX1,k,j,i) = 0.0;
                  d.dustVc[n](VX2,k,j,i) = Vk;
                  d.dustVc[n](VX3,k,j,i) = 0.0;
                }
            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}
