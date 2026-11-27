// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file test/Radiation/crossing_beams/setup.cpp
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


// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);
    real csiso = 0.1;


    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {

              d.Vc(RHO,k,j,i) = 1.e4;
              d.Vc(VX1,k,j,i) = 0.;
              d.Vc(VX2,k,j,i) = 0.;
              d.Vc(PRS,k,j,i) = d.Vc(RHO,k,j,i)*csiso*csiso;

              d.RadVc[0](ER,k,j,i) = 1.e4;
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

void InternalBoundaryRad(Fluid<RadiationPhysics> *radiation, const real t) {
  IdefixArray4D<real> Vc = radiation->Vc;
  auto *data = radiation->data;
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  int ibeg,iend,nxi,iref;

  ibeg = 0;
  iend = data->beg[IDIR];
  nxi = data->np_int[IDIR];
  iref = iend;
  idefix_for("InternalBoundaryFuncRad",
    0, data->np_tot[KDIR],
    0, data->np_tot[JDIR],
    iend, data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      //if ((j > 12) && (j < 14) && (i > 0) && (i < 4)) {

      if ((t < 100.) && (x2(j) > 0.3) && (x2(j) < 0.44) && (x1(i) > 0.5) && (x1(i) < 0.6)) {
            Vc(ER,k,j,i) = 1.e12;
            Vc(FR1,k,j,i) = 1.e12*std::cos(M_PI/4.);
            Vc(FR2,k,j,i) = 1.e12*std::sin(M_PI/4.);
      }
      if ((t < 100.) && (x2(j) > 4.56) && (x2(j) < 4.7) && (x1(i) > 0.5) && (x1(i) < 0.6)) {
            Vc(ER,k,j,i) = 1.e12;
            Vc(FR1,k,j,i) = 1.e12*std::cos(M_PI/4.);
            Vc(FR2,k,j,i) = -1.e12*std::sin(M_PI/4.);
      }

    });

}



Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output)
{
  //output.EnrollUserDefVariables(&ComputeUserVars);
  // Set the function for userdefboundary
  if(data.haveRadiation) {
    int nFrequencies = data.radiation.size();
    for(int n = 0 ; n < nFrequencies ; n++) {
      //data.radiation[n]->EnrollUserDefBoundary(&UserdefBoundaryRad);
      data.radiation[n]->EnrollInternalBoundary(&InternalBoundaryRad);
    }
    //data.hydro->EnrollUserDefBoundary(&UserdefBoundary);

  }
}
