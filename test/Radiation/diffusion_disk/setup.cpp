// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file test/Radiation/diffusion_disk/setup.cpp
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

real epsilonGlob;
real alphaGlob;
real RGlob;
real rho0Glob;
real rhominGlob;
real T0Glob;
real muGlob;
real gammaGlob;

void MySourceTerm(Hydro *hydro, const real t, const real dtin) {
  auto *data = hydro->data;
  auto units=idfx::units;

  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray4D<real> Uc = hydro->Uc;
  IdefixArray1D<real> x1=data->x[IDIR];
  IdefixArray1D<real> x2=data->x[JDIR];

  real alpha = alphaGlob;
  real R = RGlob*idfx::units.au;
  real mu = muGlob;
  real T0 = T0Glob;
  real dt=dtin;
  real Omega_K = std::sqrt(idfx::units.G*idfx::units.M_sun/(R*R*R));
  //real csiso = epsilon*R*Omega_K;
  real csiso = std::sqrt(idfx::units.k_B*T0/(mu*idfx::units.u));

  idefix_for("MySourceTerm",
    0, data->np_tot[KDIR],
    0, data->np_tot[JDIR],
    0, data->np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {

                real SE = 2.25*alpha*Omega_K*csiso*csiso*Vc(RHO,k,j,i)*units.GetDensity()/(units.GetEnergy()/units.GetTime());

                Uc(ENG,k,j,i) += dt*SE;

});
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
          Vc(ER,k,j,i) = units.ar*std::pow(T0,4.)/units.GetEnergy();
          Vc(FR1,k,j,i) = Vc(FR1,k,j,ighost);
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
          Vc(ER,k,j,i) = units.ar*std::pow(T0,4.)/units.GetEnergy();
          Vc(FR1,k,j,i) = Vc(FR1,k,j,ighost+nxi-1);
        });
    }
  }
}

void UserdefBoundary(Fluid<DefaultPhysics> *hydro, int dir, BoundarySide side, real t) {
  IdefixArray4D<real> Vc = hydro->Vc;
  auto *data = hydro->data;
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];

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
          Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost);
          Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost);
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
          Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost+nxi-1);
          Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost+nxi-1);
        });
    }
  }
}


void FluxBoundary(Fluid<DefaultPhysics> *hydro, int dir, BoundarySide side, const real t) {
    idfx::pushRegion("FluxInternal");

    auto *data = hydro->data;

    IdefixArray4D<real> FluxIDIR = data->hydro->FluxRiemann[IDIR];

    if( dir== IDIR ) {
      idefix_for("FluxInternal",
                  0, data->np_tot[KDIR],
                  0, data->np_tot[JDIR],
                  0, data->np_tot[IDIR],
         KOKKOS_LAMBDA (int k, int j, int i) {
           FluxIDIR(RHO, k, j, i) = 0.0;
           FluxIDIR(MX1, k, j, i) = 0.0;
           FluxIDIR(ENG, k, j, i) = 0.0;
      });
    }
    idfx::popRegion();

}

void FluxBoundaryOld(DataBlock & data, int dir, BoundarySide side, const real t) {
    idfx::pushRegion("FluxInternal");

    IdefixArray4D<real> FluxIDIR = data.hydro->FluxRiemann[IDIR];

    if( dir== IDIR ) {
      idefix_for("FluxInternal",
                  0, data.np_tot[KDIR],
                  0, data.np_tot[JDIR],
                  0, data.np_tot[IDIR],
         KOKKOS_LAMBDA (int k, int j, int i) {
           FluxIDIR(RHO, k, j, i) = 0.0;
           FluxIDIR(MX1, k, j, i) = 0.0;
           FluxIDIR(ENG, k, j, i) = 0.0;
      });
    }
    idfx::popRegion();

}

// Default constructor
// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output)
{
  alphaGlob = input.Get<real>("Setup","alpha",0);
  epsilonGlob = input.Get<real>("Setup","epsilon",0);
  RGlob = input.Get<real>("Setup","R",0);
  rho0Glob = input.Get<real>("Setup","rho0",0);
  rhominGlob = input.Get<real>("Setup","rhomin",0);
  T0Glob = input.Get<real>("Setup","T0",0);
  gammaGlob=input.Get<real>("Hydro","gamma",0);
  muGlob=input.Get<real>("Hydro","mu",0);

  data.hydro->EnrollUserDefBoundary(&UserdefBoundary);
  data.hydro->EnrollUserSourceTerm(&MySourceTerm);
  data.hydro->EnrollFluxBoundary(&FluxBoundary);

  // Set the function for userdefboundary
  if(data.haveRadiation) {
    int nFrequencies = data.radiation.size();
    for(int n = 0 ; n < nFrequencies ; n++) {
      data.radiation[n]->EnrollUserDefBoundary(&UserdefBoundaryRad);
    }
  }
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);

    real epsilon = epsilonGlob;
    real R = RGlob*idfx::units.au;
    real rho0 = rho0Glob;
    real rhomin = rhominGlob;
    real T0 = T0Glob;
    real mu = muGlob;
    real Omega_K = std::pow(idfx::units.G*idfx::units.M_sun/(R*R*R),0.5);
    real csiso = epsilon*R*Omega_K;
    real H2 = std::pow(epsilon*R,2.);

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {

              real x2 = std::pow(d.x[IDIR](i)*idfx::units.GetLength(),2.);

              d.Vc(RHO,k,j,i) = (rho0*std::exp(-0.5*x2/H2)+rhomin)/idfx::units.GetDensity();
              d.Vc(PRS,k,j,i) = d.Vc(RHO,k,j,i)*idfx::units.GetDensity()*idfx::units.k_B*T0/(mu*idfx::units.u)/idfx::units.GetEnergy();
              d.Vc(VX1,k,j,i) = 0.;
              d.RadVc[0](ER,k,j,i) = idfx::units.ar*std::pow(T0,4.)/idfx::units.GetEnergy();
              d.RadVc[0](FR1,k,j,i) = ZERO_F;
            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}


// Compute user variables which will be written in vtk files
void ComputeUserVars(DataBlock & data, UserDefVariablesContainer &variables) {
  // Mirror data on Host
  DataBlockHost d(data);

  // Sync it
  d.SyncFromDevice();

  // Make references to the user-defined arrays (variables is a container of IdefixHostArray3D)
  // Note that the labels should match the variable names in the input file

  for(int k = 0; k < d.np_tot[KDIR] ; k++) {
    for(int j = 0; j < d.np_tot[JDIR] ; j++) {
      for(int i = 0; i < d.np_tot[IDIR] ; i++) {

      }
    }
  }
}
