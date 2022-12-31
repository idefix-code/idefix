// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_CALCRIEMANNFLUX_HPP_
#define FLUID_CALCRIEMANNFLUX_HPP_

#if MHD == YES
#include "hlldMHD.hpp"
#include "hllMHD.hpp"
#include "roeMHD.hpp"
#include "tvdlfMHD.hpp"
#else
#include "hllcHD.hpp"
#include "hllHD.hpp"
#include "tvdlfHD.hpp"
#include "roeHD.hpp"
#endif

// Compute Riemann fluxes from states
template <typename Phys>
template <int dir>
void RiemannSolver<Phys>::CalcFlux(IdefixArray4D<real> &flux) {
  idfx::pushRegion("RiemannSolver::CalcFlux");
  if constexpr(dir == IDIR) {
    // enable shock flattening
    if(haveShockFlattening) shockFlattening.FindShock();
  }
  

  switch (mySolver) {
#if MHD == YES
    case TVDLF:
      TvdlfMHD<dir>(flux);
      break;
    case HLL:
      HllMHD<dir>(flux);
      break;
    case HLLD:
      HlldMHD<dir>(flux);
      break;
    case ROE:
      RoeMHD<dir>(flux);
      break;
#else
    case TVDLF:
      TvdlfHD<dir>(flux);
      break;
    case HLL:
      HllHD<dir>(flux);
      break;
    case HLLC:
      HllcHD<dir>(flux);
      break;
    case ROE:
      RoeHD<dir>(flux);
      break;
#endif
    default: // do nothing
      break;
  }

  idfx::popRegion();
}
#endif // FLUID_CALCRIEMANNFLUX_HPP_
