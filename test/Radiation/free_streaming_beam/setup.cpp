// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file test/Radiation/free_streaming_beam/setup.cpp
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

real csisoGlob;
real ER0Glob,ERbeamGlob,RhoGlob,AngleGlob,x1begGlob,x1endGlob,x2begGlob,x2endGlob;

void InternalBoundaryRad(Fluid<RadiationPhysics> *radiation, const real t) {
  IdefixArray4D<real> Vc = radiation->Vc;
  auto *data = radiation->data;
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  int iend;
  real x1beg = x1begGlob;
  real x1end = x1endGlob;
  real x2beg = x2begGlob;
  real x2end = x2endGlob;
  real ERbeam = ERbeamGlob;
  real Angle  = AngleGlob;

  iend = data->beg[IDIR];
  idefix_for("InternalBoundaryFuncRad",
    0, data->np_tot[KDIR],
    0, data->np_tot[JDIR],
    iend, data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      if ((x2(j) > x2beg) && (x2(j) < x2end) && (x1(i) > x1beg) && (x1(i) < x1end)) {
            Vc(ER,k,j,i) = ERbeam;;
            Vc(FR1,k,j,i) = ERbeam*std::cos(Angle);
            Vc(FR2,k,j,i) = ERbeam*std::sin(Angle);
      }
    });

}

Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output)
{
  csisoGlob = input.Get<real>("Hydro","csiso",1);
  ER0Glob = input.Get<real>("Setup","ER0",0);
  ERbeamGlob = input.Get<real>("Setup","ERbeam",0);
  AngleGlob = input.Get<real>("Setup","Angle",0);
  RhoGlob = input.Get<real>("Setup","RHO",0);
  x1begGlob = input.Get<real>("Setup","x1beg",0);
  x1endGlob = input.Get<real>("Setup","x1end",0);
  x2begGlob = input.Get<real>("Setup","x2beg",0);
  x2endGlob = input.Get<real>("Setup","x2end",0);

  if(data.haveRadiation) {
    int nFrequencies = data.radiation.size();
    for(int n = 0 ; n < nFrequencies ; n++) {
      data.radiation[n]->EnrollInternalBoundary(&InternalBoundaryRad);
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
    real csiso = csisoGlob;
    real rho0 = RhoGlob;


    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {

              d.Vc(RHO,k,j,i) = rho0;
              d.Vc(VX1,k,j,i) = ZERO_F;
              d.Vc(VX2,k,j,i) = ZERO_F;
              d.Vc(PRS,k,j,i) = d.Vc(RHO,k,j,i)*csiso*csiso;

              d.RadVc[0](ER,k,j,i) = ER0Glob;
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
