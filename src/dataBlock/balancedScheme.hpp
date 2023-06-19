// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef DATABLOCK_BALANCEDSCHEME_HPP_
#define DATABLOCK_BALANCEDSCHEME_HPP_
#include "idefix.hpp"

class DataBlock;
template<typename Phys>
class Fluid;

class BalancedScheme {
 public:
  explicit BalancedScheme(DataBlock *);
  void ComputeResidual(DataBlock *);
  void LoadResidual(DataBlock *);

  template<typename Phys>
  void SubstractBalance(Fluid<Phys>*, real);
 private:
  IdefixArray4D<real> dUc;
};

#include "fluid.hpp"
template<typename Phys>
void BalancedScheme::SubstractBalance(Fluid<Phys>* fluid, real dt) {
  idfx::pushRegion("BalancedScheme::SubstractBalance");
  auto Uc = fluid->Uc;
  auto dUc = this->dUc;

  idefix_for("delta U", 0,Uc.extent(0),
                        0,Uc.extent(1),
                        0,Uc.extent(2),
                        0,Uc.extent(3),
      KOKKOS_LAMBDA(int n, int k, int j, int i) {
        Uc(n,k,j,i) = Uc(n,k,j,i) - dUc(n,k,j,i)*dt;
      });
  idfx::popRegion();
}

#endif // DATABLOCK_BALANCEDSCHEME_HPP_
