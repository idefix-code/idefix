// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file test/Radiation/opt_thin_riemann/setup.cpp
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
real ERLGlob, ERRGlob,FR1LGlob,FR1RGlob,FR2LGlob,FR2RGlob,RhoGlob,VX1Glob,VX2Glob;

void UserdefBoundary(Fluid<DefaultPhysics> *hydro, int dir, BoundarySide side, real t) {
  IdefixArray4D<real> Vc = hydro->Vc;
  auto *data = hydro->data;

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
          Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost+nxi-1);
          Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost+nxi-1);
          Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost+nxi-1);
        });
    }
  }
}

void UserdefBoundaryRad(Fluid<RadiationPhysics> *radiation, int dir, BoundarySide side, real t) {
  IdefixArray4D<real> Vc = radiation->Vc;
  auto *data = radiation->data;

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
          Vc(FR1,k,j,i) = Vc(FR1,k,j,ighost);
          Vc(FR2,k,j,i) = Vc(FR2,k,j,ighost);
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
          Vc(ER,k,j,i) = Vc(ER,k,j,ighost+nxi-1);
          Vc(FR1,k,j,i) = Vc(FR1,k,j,ighost+nxi-1);
          Vc(FR2,k,j,i) = Vc(FR2,k,j,ighost+nxi-1);
        });
    }
  }
}

Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output)
{
  csisoGlob = input.Get<real>("Hydro","csiso",1);
  ERLGlob = input.Get<real>("Setup","ERL",0);
  ERRGlob = input.Get<real>("Setup","ERR",0);
  FR1LGlob = input.Get<real>("Setup","FR1L",0);
  FR1RGlob = input.Get<real>("Setup","FR1R",0);
  FR2LGlob = input.Get<real>("Setup","FR2L",0);
  FR2RGlob = input.Get<real>("Setup","FR2R",0);
  RhoGlob = input.Get<real>("Setup","RHO",0);
  VX1Glob = input.Get<real>("Setup","VX1",0);
  VX2Glob = input.Get<real>("Setup","VX2",0);

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

    real ERL = ERLGlob;
    real ERR= ERRGlob;
    real FR1L = FR1LGlob;
    real FR1R= FR1RGlob;
    real FR2L = FR2LGlob;
    real FR2R= FR2RGlob;
    real csiso = csisoGlob;
    real rho0 = RhoGlob;
    real vx1_0 = VX1Glob;
    real vx2_0 = VX2Glob;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {

              d.Vc(RHO,k,j,i) = rho0;
              d.Vc(VX1,k,j,i) = vx1_0;
              d.Vc(VX2,k,j,i) = vx2_0;
              d.Vc(PRS,k,j,i) = rho0*csiso*csiso;

              if (d.x[IDIR](i) < 0.){
                  d.RadVc[0](ER,k,j,i) = ERL;
                  d.RadVc[0](FR1,k,j,i) = FR1L;
                  d.RadVc[0](FR2,k,j,i) = FR2L;
              } else {
                d.RadVc[0](ER,k,j,i) = ERR;
                d.RadVc[0](FR1,k,j,i) = FR1R;
                d.RadVc[0](FR2,k,j,i) = FR2R;
              }
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
