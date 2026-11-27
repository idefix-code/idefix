// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file test/Radiation/mhd_blast/setup.cpp
//
// Last modified : 10/2026
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
real rho0Glob;
real wGlob;
real prsinGlob;
real prsoutGlob;
real B0Glob;
bool haveRadiationGlob = false;

// Default constructor
// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output)
{
  T0Glob = input.Get<real>("Setup","T0",0);
  wGlob = input.Get<real>("Setup","w",0);
  rho0Glob = input.Get<real>("Setup","rho0",0);
  prsoutGlob = input.Get<real>("Setup","prs_out",0);
  prsinGlob = input.Get<real>("Setup","prs_in",0);
  B0Glob = input.Get<real>("Setup","B0",0);
  if (input.CheckBlock("Rad")) {
    haveRadiationGlob = true;
  }
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);
    real T0 = T0Glob;
    real rho0 = rho0Glob;
    real w = wGlob;
    real r2,T;
    real prs_out = prsoutGlob;
    real prs_in = prsinGlob;
    real B0 = B0Glob;

    real unit_density = idfx::units.GetDensity();
    real unit_velocity = idfx::units.GetVelocity();
    real unit_energy = idfx::units.GetEnergy();

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
              if (GEOMETRY==CARTESIAN){
                r2 = d.x[IDIR](i)*d.x[IDIR](i)+d.x[JDIR](j)*d.x[JDIR](j)+d.x[KDIR](k)*d.x[KDIR](k);
              } else if (GEOMETRY==SPHERICAL) {
                r2 = d.x[IDIR](i)*d.x[IDIR](i);
              }
              d.Vc(RHO,k,j,i) = rho0;
              d.Vc(PRS,k,j,i) = prs_out;
              if (r2 < w*w) d.Vc(PRS,k,j,i) = prs_in;
              d.Vc(VX1,k,j,i) = ZERO_F;
              d.Vc(VX2,k,j,i) = ZERO_F;
              d.Vc(VX3,k,j,i) = ZERO_F;
              T = prs_in/rho0*idfx::units.GetKelvin();
              if (haveRadiationGlob) {
                d.RadVc[0](ER,k,j,i) = idfx::units.ar*std::pow(T,4)/unit_energy;
                d.RadVc[0](FR1,k,j,i) = ZERO_F;
                d.RadVc[0](FR2,k,j,i) = ZERO_F;
                d.RadVc[0](FR3,k,j,i) = ZERO_F;
              }
              #if MHD == YES
              d.Vs(BX1s,k,j,i) = B0;
              d.Vs(BX2s,k,j,i) = 0.0;
              d.Vc(BX3,k,j,i) = 0.0;
              #endif
            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}
