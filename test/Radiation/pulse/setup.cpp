// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file test/Radiation/pulse/setup.cpp
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


real T0Glob;
real muGlob;
real rho0Glob;
real wGlob;

// Default constructor
// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output)
{
  T0Glob = input.Get<real>("Setup","T0",0);
  wGlob = input.Get<real>("Setup","w",0);
  rho0Glob = input.Get<real>("Setup","rho0",0);
  muGlob = input.Get<real>("Hydro","mu",0);
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);
    real mu = muGlob;
    real T0 = T0Glob;
    real rho0 = rho0Glob;
    real w = wGlob;
    real r2,T;


    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
              if (GEOMETRY==CARTESIAN){
                r2 = EXPAND(d.x[IDIR](i)*d.x[IDIR](i),+d.x[JDIR](j)*d.x[JDIR](j),+d.x[KDIR](k)*d.x[KDIR](k));
              } else if (GEOMETRY==SPHERICAL) {
                r2 = d.x[IDIR](i)*d.x[IDIR](i);
              }
              d.Vc(RHO,k,j,i) = rho0;
              d.Vc(PRS,k,j,i) = d.Vc(RHO,k,j,i)*T0/idfx::units.GetKelvin();
              EXPAND( d.Vc(VX1,k,j,i) = ZERO_F;,
                      d.Vc(VX2,k,j,i) = ZERO_F;,
                      d.Vc(VX3,k,j,i) = ZERO_F;)
              T = T0*(1.+100.*std::exp(-r2/(w*w)));
              d.RadVc[0](ER,k,j,i) = idfx::units.ar*std::pow(T,4)/idfx::units.GetEnergy();
              EXPAND(d.RadVc[0](FR1,k,j,i) = ZERO_F;,
                     d.RadVc[0](FR2,k,j,i) = ZERO_F;,
                     d.RadVc[0](FR3,k,j,i) = ZERO_F;)
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
