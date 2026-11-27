// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file test/Radiation/novae_like/setup.cpp
//
// Last modified : 02/2026
//
// Copyright(C) by :
// - Nicolas Scepi <nicolas.scepi@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "idefix.hpp"
#include "setup.hpp"

real R0Glob;
real rho0Glob;
real rhominGlob;
real T0Glob;
real muGlob;
real gammaGlob;
real alphaMRIGlob;
real MGlob;
real GGlob;
real densityFloorGlob;
real rhoindexGlob;
real x1begGlob;


void MyViscosity(DataBlock &data, const real t, IdefixArray3D<real> &eta1, IdefixArray3D<real> &eta2) {
  auto units = idfx::units;

  IdefixArray3D<real> InvDt = data.hydro->InvDt;
  IdefixArray4D<real> Vc=data.hydro->Vc;
  IdefixArray1D<real> r=data.x[IDIR];
  IdefixArray1D<real> th=data.x[JDIR];
  real alphaMRI = alphaMRIGlob;
  real CG = MGlob*GGlob;
  real R0 = R0Glob;

  idefix_for("MyViscosity",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real R = FMAX(r(i)*sin(th(j)),R0);
                real cs2 = Vc(PRS,k,j,i)/Vc(RHO,k,j,i);
                real Omega = std::sqrt(CG)*std::pow(R,-1.5);
                real alpha = alphaMRI;
                eta1(k,j,i) = alpha*cs2*Vc(RHO,k,j,i)/Omega;
                eta2(k,j,i) = 0.;
              });

}


// Compute user variables which will be written in vtk files
void ComputeUserVars(DataBlock & data, UserDefVariablesContainer &variables) {
  // Mirror data on Host
  DataBlockHost d(data);

  // Sync it
  d.SyncFromDevice();

  auto units = idfx::units;

  IdefixHostArray1D<real> x1=d.x[IDIR];
  IdefixHostArray1D<real> x1l=d.xl[IDIR];
  IdefixHostArray1D<real> x2=d.x[JDIR];
  IdefixHostArray1D<real> x2l=d.xl[JDIR];
  IdefixHostArray3D<real> A1=d.A[IDIR];
  IdefixHostArray3D<real> A2=d.A[JDIR];
  IdefixHostArray3D<real> dV=d.dV;

  IdefixArray4D<real> Vc=(&data)->hydro->Vc;

  IdefixHostArray3D<real> divF  = variables["divF"];
  IdefixHostArray3D<real> A1_out  = variables["A1"];

  Kokkos::deep_copy(variables["dV"], dV);
  Kokkos::deep_copy(variables["divF"], divF);
  Kokkos::deep_copy(variables["A1"], A1_out);

}


void UserdefBoundaryNoStress(Fluid<DefaultPhysics> *hydro, int dir, BoundarySide side, real t) {
  IdefixArray4D<real> Vc = hydro->Vc;
  auto *data = hydro->data;
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  real rhomin = densityFloorGlob/idfx::units.GetDensity();

  if(dir==IDIR) {
    int ighost,nxi,iend,ibeg;
    if(side == left) {
      ighost = data->nghost[IDIR];
      ibeg = 0;
      iend = data->beg[IDIR];
      idefix_for("UserDefBoundary",
        0, data->np_tot[KDIR],
        0, data->np_tot[JDIR],
        ibeg, iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost);
          Vc(VX3,k,j,i) = x1(i)*std::sin(x2(j))*Vc(VX3,k,j,ighost)/(x1(ighost)*std::sin(x2(j)));

          Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost)/Vc(RHO,k,j,ighost)*Vc(RHO,k,j,i);
          if(Vc(VX1,k,j,ighost)>=ZERO_F){
            Vc(VX1,k,j,i)=ZERO_F;
          }else {
            Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost);
          }
          Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost);
        });
    } else if (side ==right) {
      ighost = data->nghost[IDIR];
      nxi = data->np_int[IDIR];
      ibeg = data->end[IDIR];
      iend = data->np_tot[IDIR];
      idefix_for("UserDefBoundary",
        0, data->np_tot[KDIR],
        0, data->np_tot[JDIR],
        ibeg, iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost+nxi-1);
          Vc(VX3,k,j,i) = Vc(VX3,k,j,ighost+nxi-1);

          Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost+nxi-1)/Vc(RHO,k,j,ighost+nxi-1)*Vc(RHO,k,j,i);
          if(Vc(VX1,k,j,ighost+nxi-1)<=ZERO_F){
            Vc(VX1,k,j,i)=ZERO_F;
          }else {
            Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost+nxi-1);
          }
          Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost+nxi-1);
        });
    }
  }
}

void UserdefBoundaryRad(Fluid<RadiationPhysics> *radiation, int dir, BoundarySide side, real t) {
  IdefixArray4D<real> Vc = radiation->Vc;
  auto *data = radiation->data;
  auto units=idfx::units;

  real T0 = T0Glob;

  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  if(dir==IDIR) {
    int ighost,nxi,iend,ibeg;
    if(side == left) {
      ighost = data->nghost[IDIR];
      ibeg = 0;
      iend = data->beg[IDIR];
      idefix_for("UserDefBoundaryRad",
        0, data->np_tot[KDIR],
        0, data->np_tot[JDIR],
        ibeg, iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vc(ER,k,j,i) = Vc(ER,k,j,ighost);
          if (Vc(FR1,k,j,ighost) >=ZERO_F){
            Vc(FR1,k,j,i) = ZERO_F;
          } else {
            Vc(FR1,k,j,i) = Vc(FR1,k,j,ighost);
          }
          Vc(FR2,k,j,i) = Vc(FR2,k,j,ighost);
          Vc(FR3,k,j,i) = Vc(FR3,k,j,ighost);
        });
    } else if (side==right){
      ighost = data->nghost[IDIR];
      nxi = data->np_int[IDIR];
      ibeg = data->end[IDIR];
      iend =data->np_tot[IDIR];
      idefix_for("UserDefBoundaryRad",
        0, data->np_tot[KDIR],
        0, data->np_tot[JDIR],
        ibeg, iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vc(ER,k,j,i) = units.ar*std::pow(10.,4.)/units.GetEnergy();
          if (Vc(FR1,k,j,ighost+nxi-1) <=ZERO_F){
            Vc(FR1,k,j,i) = ZERO_F;
          } else {
            Vc(FR1,k,j,i) = Vc(FR1,k,j,ighost+nxi-1);
          }
          Vc(FR2,k,j,i) = Vc(FR2,k,j,ighost+nxi-1);
          Vc(FR3,k,j,i) = Vc(FR3,k,j,ighost+nxi-1);
        });
    }
  }
}


void InternalBoundary(Hydro *hydro, const real t) {
  auto *data = hydro->data;
  auto units = idfx::units;

  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray4D<real> Vs = hydro->Vs;
  IdefixArray1D<real> x1=data->x[IDIR];
  IdefixArray1D<real> x2=data->x[JDIR];

  real densityFloor = densityFloorGlob;
  real rhoindex = rhoindexGlob;
  real x1beg = x1begGlob;

  idefix_for("InternalBoundary",
    0, data->np_tot[KDIR],
    0, data->np_tot[JDIR],
    0, data->np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real densfloor = densityFloor*std::pow(x1beg/x1(i),rhoindex);
                if(Vc(RHO,k,j,i)*units.GetDensity() <= densfloor) {
                  real T= Vc(PRS,k,j,i)/Vc(RHO,k,j,i);
                  Vc(RHO,k,j,i)=densfloor/units.GetDensity();
                  Vc(PRS,k,j,i)=T*Vc(RHO,k,j,i);
                  Vc(VX1,k,j,i)=ZERO_F;
                  Vc(VX2,k,j,i)=ZERO_F;
                }
              });
}

// Default constructor
// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output)
{
  auto Vc = data.hydro->Vc;

  // Mirror data on Host
  alphaMRIGlob = input.Get<real>("Setup","alphaMRI",0);
  R0Glob = input.Get<real>("Setup","R0",0);
  rho0Glob = input.Get<real>("Setup","rho0",0);
  rhoindexGlob = input.Get<real>("Setup","rhoindex",0);
  rhominGlob = input.Get<real>("Setup","rhomin",0);
  T0Glob = input.Get<real>("Setup","T0",0);
  densityFloorGlob = input.Get<real>("Setup","density_floor",0);

  GGlob = input.Get<real>("Gravity","gravCst",0);
  MGlob = input.Get<real>("Gravity","Mcentral",0);

  muGlob = input.Get<real>("Hydro","mu",0);
  gammaGlob=data.hydro->eos->GetGamma();

  x1begGlob = grid.xbeg[IDIR];

  data.hydro->EnrollInternalBoundary(&InternalBoundary);
  data.hydro->viscosity->EnrollViscousDiffusivity(&MyViscosity);
  data.hydro->EnrollUserDefBoundary(&UserdefBoundaryNoStress);
  data.radiation[0]->EnrollUserDefBoundary(&UserdefBoundaryRad);

}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);

    auto units=idfx::units;

    real R0 = R0Glob;
    real rho0 = rho0Glob;
    real densityFloor = densityFloorGlob;
    real rhoindex = rhoindexGlob;
    real T0 = T0Glob;
    real mu = muGlob;
    real CG = MGlob*GGlob;
    real x1beg = x1begGlob;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {

              real R = FMAX(d.x[IDIR](i)*std::sin(d.x[JDIR](j)),R0);
              real z2 = std::pow(d.x[IDIR](i)*std::cos(d.x[JDIR](j)),2.);
              real Omega = std::sqrt(CG)*std::pow(R,-1.5);
              real cs = std::sqrt(T0/units.GetKelvin()/mu);
              real H = cs/Omega;

              real densfloor = densityFloor*std::pow(x1beg/d.x[IDIR](i),rhoindex);
              d.Vc(RHO,k,j,i) = (rho0*std::pow(R0/R,rhoindex)*std::exp(-0.25*M_PI*z2/(H*H))+densfloor)/idfx::units.GetDensity();
              d.Vc(PRS,k,j,i) = d.Vc(RHO,k,j,i)*T0/units.GetKelvin()/mu;
              d.Vc(VX1,k,j,i) = 0.;
              d.Vc(VX2,k,j,i) = 0.;
              d.Vc(VX3,k,j,i) =  Omega*R;

              real T = d.Vc(PRS,k,j,i)/d.Vc(RHO,k,j,i)*units.GetKelvin()*mu;
              d.RadVc[0](ER,k,j,i) = idfx::units.ar*std::pow(T,4.)/idfx::units.GetEnergy();
              d.RadVc[0](FR1,k,j,i) = ZERO_F;
              d.RadVc[0](FR2,k,j,i) = ZERO_F;
              d.RadVc[0](FR3,k,j,i) = ZERO_F;

            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}
