// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file test/Radiation/shocks/setup.cpp
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

real rho0Glob,vx1Glob,T0Glob,muGlob;

void UserdefBoundary(Fluid<DefaultPhysics> *hydro, int dir, BoundarySide side, real t) {
  IdefixArray4D<real> Vc = hydro->Vc;
  auto *data = hydro->data;

  if(dir==IDIR) {
    int ighost,nxi,iend,ibeg;
    if (side == right) {
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

void UserdefBoundaryRad(Fluid<RadiationPhysics> *radiation, int dir, BoundarySide side, real t) {
  IdefixArray4D<real> Vc = radiation->Vc;
  auto *data = radiation->data;

  if(dir==IDIR) {
    int ighost,nxi,iend,ibeg;
    if (side==right){
      ighost = data->nghost[IDIR];
      nxi = data->np_int[IDIR];
      ibeg = data->end[IDIR];
      iend =data->np_tot[IDIR];
      idefix_for("UserDefBoundaryRad",
        0, data->np_tot[KDIR],
        0, data->np_tot[JDIR],
        ibeg, iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vc(ER,k,j,i) = Vc(ER,k,j,ighost+nxi-1);
          Vc(FR1,k,j,i) = Vc(FR1,k,j,ighost+nxi-1);
        });
    }
  }
}

// Default constructor
// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output)
{
  rho0Glob = input.Get<real>("Setup","rho0",0);
  vx1Glob = input.Get<real>("Setup","vx1",0);
  T0Glob = input.Get<real>("Setup","T0",0);
  muGlob = input.Get<real>("Hydro","mu",0);

  // Set the function for userdefboundary
  data.hydro->EnrollUserDefBoundary(&UserdefBoundary);
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

    real T = T0Glob;
    real rho0 = rho0Glob;
    real vx1 = vx1Glob;
    real mu = muGlob;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
              d.Vc(RHO,k,j,i) = rho0/idfx::units.GetDensity();
              d.Vc(PRS,k,j,i) = d.Vc(RHO,k,j,i)*idfx::units.GetDensity()*idfx::units.k_B*T/(mu*idfx::units.u*idfx::units.GetEnergy());
              d.Vc(VX1,k,j,i) = vx1/idfx::units.GetVelocity();
              d.RadVc[0](ER,k,j,i) = idfx::units.ar*std::pow(T,4)/(idfx::units.GetEnergy());
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
