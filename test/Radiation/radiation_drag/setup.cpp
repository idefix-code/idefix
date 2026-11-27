// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file test/Radiation/radiation_drag/setup.cpp
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

real BigCGlob;
real BigPGlob;
real T0Glob;
real muGlob;
real gammaGlob;



void FluxBoundary(Fluid<DefaultPhysics> *hydro, int dir, BoundarySide side, const real t) {
    idfx::pushRegion("FluxInternal");

    auto *data = hydro->data;

    IdefixArray4D<real> FluxIDIR = data->hydro->FluxRiemann[IDIR];

    idefix_for("FluxInternal",
                0, data->np_tot[KDIR],
                0, data->np_tot[JDIR],
                0, data->np_tot[IDIR],
       KOKKOS_LAMBDA (int k, int j, int i) {
         FluxIDIR(RHO, k, j, i) = 0.0;
    });
    idfx::popRegion();

}


void FluxBoundaryRad(Fluid<RadiationPhysics> *radiation, int dir, BoundarySide side, const real t) {
    idfx::pushRegion("FluxInternal");

    auto *data = radiation->data;

    IdefixArray4D<real> FluxIDIR = data->radiation[0]->FluxRiemann[IDIR];

    idefix_for("FluxInternal",
               0, data->np_tot[KDIR],
               0, data->np_tot[JDIR],
               0, data->np_tot[IDIR],
      KOKKOS_LAMBDA (int k, int j, int i) {
        FluxIDIR(ER, k, j, i) = 0.0;
    });

    idfx::popRegion();

}

void InternalBoundary(Fluid<DefaultPhysics> *hydro, const real t) {
  IdefixArray4D<real> Vc = hydro->Vc;
  auto *data = hydro->data;

  real mu = muGlob;
  real rho = idfx::units.ar/BigPGlob*std::pow(mu*idfx::units.u/idfx::units.k_B,4.)*std::pow(idfx::units.c/BigCGlob,6.)/idfx::units.GetDensity();

  idefix_for("InternalBoundary",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
        Vc(RHO,k,j,i) = rho;
//        std::printf("rho=%e in InternalBoundary\n",Vc(RHO,k,j,i));
    });

}

void InternalBoundaryRad(Fluid<RadiationPhysics> *radiation, const real t) {
  IdefixArray4D<real> Vc = radiation->Vc;
  auto *data = radiation->data;

  real mu = muGlob;
  real cs = idfx::units.c/BigCGlob/idfx::units.GetVelocity();
  real rho = idfx::units.ar/BigPGlob*std::pow(mu*idfx::units.u/idfx::units.k_B,4.)*std::pow(idfx::units.c/BigCGlob,6.)/idfx::units.GetDensity();
  real P = rho*cs*cs;
  real Er = BigPGlob*P;

  idefix_for("InternalBoundary",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
        Vc(ER,k,j,i) = Er;
  //      std::printf("Er=%e in InternalBoundary\n",Vc(ER,k,j,i));
    });

}

// Default constructor
// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output)
{
  gammaGlob=input.Get<real>("Hydro","gamma",0);
  muGlob=input.Get<real>("Hydro","mu",0);
  BigCGlob=input.Get<real>("Setup","BigC",0);
  BigPGlob=input.Get<real>("Setup","BigP",0);

  data.hydro->EnrollFluxBoundary(&FluxBoundary);
  data.hydro->EnrollInternalBoundary(&InternalBoundary);

  // Set the function for userdefboundary
  if(data.haveRadiation) {
    int nFrequencies = data.radiation.size();
    for(int n = 0 ; n < nFrequencies ; n++) {
      data.radiation[n]->EnrollFluxBoundary(&FluxBoundaryRad);
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

    real mu = muGlob;
    real cs = idfx::units.c/BigCGlob/idfx::units.GetVelocity();
    real rho = idfx::units.ar/BigPGlob*std::pow(mu*idfx::units.u/idfx::units.k_B,4.)*std::pow(idfx::units.c/BigCGlob,6.)/idfx::units.GetDensity();
    real P = rho*cs*cs;
    real Er = BigPGlob*P;
    real v = 9.480e9/idfx::units.GetVelocity();
    std::printf("rho=%e P=%e cs=%e Er=%e v=%e\n",rho,P,cs,Er,v);


    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {

              d.Vc(PRS,k,j,i) = P;
              d.Vc(RHO,k,j,i) = rho;
              d.Vc(VX1,k,j,i) = v;
              d.RadVc[0](ER,k,j,i) = Er;
              //d.RadVc[0](ER,k,j,i) = ZERO_F;
              d.RadVc[0](FR1,k,j,i) = ZERO_F;
            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}
