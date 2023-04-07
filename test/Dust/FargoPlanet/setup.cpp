#include <algorithm>
#include "idefix.hpp"
#include "setup.hpp"

real sigmaSlopeGlob;
real h0Glob;
real betaGlob;
real HidealGlob;
real AmMidGlob;
real gammaGlob;
real densityFloorGlob;
real alphaGlob;
std::vector<real> taus;


void MySoundSpeed(DataBlock &data, const real t, IdefixArray3D<real> &cs) {
  IdefixArray1D<real> x1=data.x[IDIR];
  real h0 = h0Glob;
  idefix_for("MySoundSpeed",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real R = x1(i);
                cs(k,j,i) = h0/sqrt(R);
              });
}

void DustDrag(Fluid<DustPhysics> *dust, const real t, const real dt) {
  IdefixArray4D<real> dustVc = dust->Vc;
  IdefixArray4D<real> dustUc = dust->Uc;
  auto data = dust->data;
  IdefixArray4D<real> Vc = data->hydro->Vc;
  IdefixArray4D<real> Uc = data->hydro->Uc;
  IdefixArray1D<real> x1=data->x[IDIR];

  real tau = taus[dust->instanceNumber];
  if(tau < 2*dt) tau=2*dt;
  real h0 = h0Glob;
  idefix_for("DustDrag",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                // Dust drag
                real cs = h0/sqrt(x1(i));

                real tauloc = tau*h0/(cs*Vc(RHO,k,j,i)); // So that stopping time=tau at R=1 and rho=1

                for(int n = MX1 ; n < MX3 ; n++) {
                  real dp = dt*dustVc(RHO,k,j,i) * (dustVc(n,k,j,i) - Vc(n,k,j,i)) / tau;
                  Uc(n,k,j,i) += dp;
                  dustUc(n,k,j,i) -= dp;
                }
              });
  }


// User-defined boundaries
void UserdefBoundary(DataBlock& data, int dir, BoundarySide side, real t) {
  IdefixArray4D<real> Vc = data.hydro->Vc;
  IdefixArray1D<real> x1 = data.x[IDIR];
  IdefixArray1D<real> x3 = data.x[KDIR];
  if(dir==IDIR) {
    int ighost,ibeg,iend;
    if(side == left) {
      ighost = data.beg[IDIR];
      ibeg = 0;
      iend = data.beg[IDIR];
      idefix_for("UserDefBoundary",0,data.np_tot[KDIR],0,data.np_tot[JDIR],ibeg,iend,
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
      ighost = data.end[IDIR]-1;
      ibeg=data.end[IDIR];
      iend=data.np_tot[IDIR];
      idefix_for("UserDefBoundary",0,data.np_tot[KDIR],0,data.np_tot[JDIR],ibeg,iend,
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
    taus = std::vector<real>(nSpecies);
    for(int n = 0 ; n < nSpecies ; n++) {
      data.dust[n]->EnrollUserDefBoundary(&UserdefBoundaryDust);
      data.dust[n]->EnrollUserSourceTerm(&DustDrag);
      taus[n] = input.Get<real>("Setup","taus",n);
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
