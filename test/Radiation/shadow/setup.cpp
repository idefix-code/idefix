// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file test/Radiation/shadow/setup.cpp
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

real rho0Glob;
real rho1Glob;
real T0Glob;
real TinjGlob;
real muGlob;
real x02Glob;
real y02Glob;


void UserdefBoundaryRad(Fluid<RadiationPhysics> *radiation, int dir, BoundarySide side, real t) {
  IdefixArray4D<real> Vc = radiation->Vc;
  auto *data = radiation->data;
  auto units=idfx::units;

  real Tinj = TinjGlob;

  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  if(dir==IDIR) {
    int ibeg,iend;
    if(side == left) {
      ibeg = 0;
      iend = data->beg[IDIR];
      idefix_for("UserDefBoundaryRad",
        0, data->np_tot[KDIR],
        0, data->np_tot[JDIR],
        ibeg, iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vc(ER,k,j,i) = units.ar*std::pow(Tinj,4)/units.GetEnergy();
          Vc(FR1,k,j,i) = Vc(ER,k,j,i);
          Vc(FR2,k,j,i) = ZERO_F;
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
    int ighost,ibeg,iend;
    if(side == left) {
      ighost = data->beg[IDIR];
      ibeg = 0;
      iend = data->beg[IDIR];
      idefix_for("UserDefBoundary",
        0, data->np_tot[KDIR],
        0, data->np_tot[JDIR],
        ibeg, iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vc(RHO,k,j,i) = Vc(RHO,k,j,2*ighost-i-1);
          Vc(PRS,k,j,i) = Vc(PRS,k,j,2*ighost-i-1);
          Vc(VX1,k,j,i) = -Vc(VX1,k,j,2*ighost-i-1);
          Vc(VX2,k,j,i) = -Vc(VX2,k,j,2*ighost-i-1);
        });
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

    real C_ar = idfx::units.ar;
    real mu = muGlob;
    real rho0 = rho0Glob;
    real rho1 = rho1Glob;
    real T0 = T0Glob;
    real x02 = x02Glob/(idfx::units.GetLength()*idfx::units.GetLength());
    real y02 = y02Glob/(idfx::units.GetLength()*idfx::units.GetLength());

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
              real delta = 10.*(d.x[IDIR](i)*d.x[IDIR](i)/x02+d.x[JDIR](j)*d.x[JDIR](j)/y02-1.);
              d.Vc(RHO,k,j,i) = rho0+(rho1-rho0)/(1.+std::exp(delta));
              d.Vc(PRS,k,j,i) = d.Vc(RHO,k,j,i)*T0/idfx::units.GetKelvin();
              d.Vc(VX1,k,j,i) = ZERO_F;
              d.Vc(VX2,k,j,i) = ZERO_F;
              d.RadVc[0](ER,k,j,i) = C_ar*std::pow(T0,4)/idfx::units.GetEnergy();
              d.RadVc[0](FR1,k,j,i) = ZERO_F;
              d.RadVc[0](FR2,k,j,i) = ZERO_F;
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

// Default constructor
// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output)
{
  muGlob = input.Get<real>("Hydro","mu",0);
  rho0Glob = input.Get<real>("Rad","kappa",3);
  T0Glob = input.Get<real>("Rad","kappa",4);
  rho1Glob = input.Get<real>("Setup","rho1",0);
  TinjGlob = input.Get<real>("Setup","Tinj",0);
  x02Glob = input.Get<real>("Setup","x02",0);
  y02Glob = input.Get<real>("Setup","y02",0);
  //output.EnrollUserDefVariables(&ComputeUserVars);
  // Set the function for userdefboundary
  if(data.haveRadiation) {
    int nFrequencies = data.radiation.size();
    for(int n = 0 ; n < nFrequencies ; n++) {
      data.radiation[n]->EnrollUserDefBoundary(&UserdefBoundaryRad);
    }
    data.hydro->EnrollUserDefBoundary(&UserdefBoundary);
    //data.hydro->EnrollInternalBoundary(&InternalBoundary);

  }
}
