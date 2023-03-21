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
real tauGlob;


void MySoundSpeed(DataBlock &data, const real t, IdefixArray3D<real> &cs) {
  IdefixArray1D<real> x1=data.x[IDIR];
  real h0 = h0Glob;
  idefix_for("MySoundSpeed",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real R = x1(i);
                cs(k,j,i) = h0/sqrt(R);
              });
}

void DustFeedback(DataBlock &data, const real t, const real dt) {
  real tau0 = tauGlob;
  auto Uc = data.hydro->Uc;
  auto Vc = data.hydro->Vc;
  auto dustUc = data.dust[0]->Uc;
  auto dustVc = data.dust[0]->Vc;
  idefix_for("MySourceTerm",MX1,MX1+COMPONENTS,0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int n, int k, int j, int i) {
                real tau = tau0 / Vc(RHO,k,j,i);
                // Avoid too small stopping times
                if(tau<2*dt) tau=2*dt;

                real dp = dt*dustVc(RHO,k,j,i) * (dustVc(n,k,j,i) - Vc(n,k,j,i)) / tau;
                Uc(n,k,j,i) += dp;
                dustUc(n,k,j,i) -= dp;

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

void UserdefDustBoundary(DataBlock& data, int dir, BoundarySide side, real t) {
  IdefixArray4D<real> Vc = data.dust[0]->Vc;
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
  data.dust[0]->EnrollUserDefBoundary(&UserdefDustBoundary);
  data.hydro->EnrollIsoSoundSpeed(&MySoundSpeed);
  data.hydro->EnrollUserSourceTerm(&DustFeedback);
  if(data.haveFargo)
    data.fargo->EnrollVelocity(&FargoVelocity);
  sigmaSlopeGlob = input.Get<real>("Setup","sigmaSlope",0);
  h0Glob = input.Get<real>("Setup","h0",0);
  tauGlob = input.Get<real>("Setup","taus",0); // stopping time @ R=1
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

                d.dustVc[0](RHO,k,j,i) = 1e-2*d.Vc(RHO,k,j,i);
                d.dustVc[0](VX1,k,j,i) = 0.0;
                d.dustVc[0](VX2,k,j,i) = Vk;
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
