// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>
#include "../idefix.hpp"
#include "dataBlock.hpp"
#include "dump.hpp"
#include "balancedScheme.hpp"

BalancedScheme::BalancedScheme(DataBlock &data) {
  if(data.haveDust) {
    IDEFIX_WARNING("The balanced scheme does not apply to dust fields");
  }
  auto Uc = data.hydro->Uc;
  IdefixArray4D<real> dUc("Balance::dUc",
                          Uc.extent(0),
                          Uc.extent(1),
                          Uc.extent(2),
                          Uc.extent(3));
}

void BalancedScheme::ComputeResidual(DataBlock &data) {
  data.ResetStage();
  data.SetBoundaries();
  if(data.haveFargo) data.fargo->SubstractVelocity(data.t);
  data.PrimToCons();
  // Backup initial state
  Kokkos::deep_copy(UcInit,data.hydro->Uc):

  // Evolve the stage by one timeStep
  if(data.haveGravity) {
    data.gravity->ComputeGravity(ncycles);
  }
  data.EvolveStage();

  // Compute the variation of conserved quantities
  real dt = data.dt;

  idefix_for("delta U", 0,Uc.extent(0),
                        0,Uc.extent(1),
                        0,Uc.extent(2),
                        0,Uc.extent(3),
      KOKKOS_LAMBDA(int n, int k, int j, int i) {
        dUc(n,k,j,i) = (Uc(n,k,j,i) - dUc(n,k,j,i))/dt;
      });

  Dump dump(&data);

  // Register all of the variables
  for(int i = 0 ; i < dUc.extent(0) ;  i++) {
    dump->RegisterVariable(dUc,std::string("Res")+std::to_string(i),i);
  }
  dump.SetFilename("residual");
  dump.write()
}


void BalancedScheme::LoadResidual(DataBlock &data) {
  if(data.haveDust) {
    IDEFIX_WARNING("The balanced scheme does not apply to dust fields");
  }

  Dump dump(&data);

  // Register all of the variables
  for(int i = 0 ; i < dUc.extent(0) ;  i++) {
    dump->RegisterVariable(dUc,std::string("Res")+std::to_string(i),i);
  }
  dump.SetFilename("residual");
  dump.read()
}

IdefixArray4D<real> BalancedScheme::GetResidual() {
  return this->dUc;
}
