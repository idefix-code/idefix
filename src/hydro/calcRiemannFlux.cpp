// ***********************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************************

#include "../idefix.hpp"
#include "hydro.hpp"

#if MHD == YES
#include "solversMHD.hpp"
#else
#include "solversHD.hpp"
#endif

// Compute Riemann fluxes from states
void Hydro::CalcRiemannFlux(int dir, const real t) {
  idfx::pushRegion("Hydro::CalcRiemannFlux");

  if(this->haveHall == UserDefFunction && dir == IDIR) {
    if(this->hallDiffusivityFunc)
      hallDiffusivityFunc(*data, t, this->xHall);
    else
      IDEFIX_ERROR("No user-defined Hall diffusivity function has been enrolled");
  }

  switch (mySolver) {
#if MHD == YES
    case TVDLF:
      switch(dir) {
        case IDIR:
          this->TvdlfMHD<IDIR,ARG_EXPAND(MX1,MX2,MX3),
                   ARG_EXPAND(BX1,BX2,BX3)>();
          break;
  #if DIMENSIONS >= 2
        case JDIR:
          this->TvdlfMHD<JDIR,ARG_EXPAND(MX2,MX1,MX3),
                   ARG_EXPAND(BX2,BX1,BX3)>();
          break;
  #endif
  #if DIMENSIONS == 3
        case KDIR:
          this->TvdlfMHD<KDIR,ARG_EXPAND(MX3,MX1,MX2),
                   ARG_EXPAND(BX3,BX1,BX2)>();
          break;
  #endif
      }
      break;
    case HLL:
      switch(dir) {
        case IDIR:
          this->HllMHD<IDIR,ARG_EXPAND(MX1,MX2,MX3),
                 ARG_EXPAND(BX1,BX2,BX3)>();
          break;
  #if DIMENSIONS >= 2
        case JDIR:
          this->HllMHD<JDIR,ARG_EXPAND(MX2,MX1,MX3),
                 ARG_EXPAND(BX2,BX1,BX3)>();
          break;
  #endif
  #if DIMENSIONS == 3
        case KDIR:
          this->HllMHD<KDIR,ARG_EXPAND(MX3,MX1,MX2),
                 ARG_EXPAND(BX3,BX1,BX2)>();
          break;
  #endif
      }
      break;
    case HLLD:
      switch(dir) {
        case IDIR:
          this->HlldMHD<IDIR,ARG_EXPAND(MX1,MX2,MX3),
                  ARG_EXPAND(BX1,BX2,BX3)>();
          break;
  #if DIMENSIONS >= 2
        case JDIR:
          this->HlldMHD<JDIR,ARG_EXPAND(MX2,MX1,MX3),
                  ARG_EXPAND(BX2,BX1,BX3)>();
          break;
  #endif
  #if DIMENSIONS == 3
        case KDIR:
          this->HlldMHD<KDIR,ARG_EXPAND(MX3,MX1,MX2),
                  ARG_EXPAND(BX3,BX1,BX2)>();
          break;
  #endif
      }
      break;
    case ROE:
      switch(dir) {
        case IDIR:
          this->RoeMHD<IDIR,ARG_EXPAND(MX1,MX2,MX3),
                 ARG_EXPAND(BX1,BX2,BX3)>();
          break;
  #if DIMENSIONS >= 2
        case JDIR:
          this->RoeMHD<JDIR,ARG_EXPAND(MX2,MX1,MX3),
                 ARG_EXPAND(BX2,BX1,BX3)>();
          break;
  #endif
  #if DIMENSIONS == 3
        case KDIR:
          this->RoeMHD<KDIR,ARG_EXPAND(MX3,MX1,MX2),
                 ARG_EXPAND(BX3,BX1,BX2)>();
          break;
  #endif
      }
      break;
#else
    case TVDLF:
      switch(dir) {
        case IDIR:
          this->TvdlfHD<IDIR,MX1,MX2,MX3>();
          break;
  #if DIMENSIONS >= 2
        case JDIR:
          this->TvdlfHD<JDIR,MX2,MX1,MX3>();
          break;
  #endif
  #if DIMENSIONS == 3
        case KDIR:
          this->TvdlfHD<KDIR,MX3,MX1,MX2>();
          break;
  #endif
      }
      break;
    case HLL:
      switch(dir) {
        case IDIR:
          this->HllHD<IDIR,MX1,MX2,MX3>();
          break;
  #if DIMENSIONS >= 2
        case JDIR:
          this->HllHD<JDIR,MX2,MX1,MX3>();
          break;
  #endif
  #if DIMENSIONS == 3
        case KDIR:
          this->HllHD<KDIR,MX3,MX1,MX2>();
          break;
  #endif
      }
      break;
    case HLLC:
      switch(dir) {
        case IDIR:
          this->HllcHD<IDIR,MX1,MX2,MX3>();
          break;
  #if DIMENSIONS >= 2
        case JDIR:
          this->HllcHD<JDIR,MX2,MX1,MX3>();
          break;
  #endif
  #if DIMENSIONS == 3
        case KDIR:
          this->HllcHD<KDIR,MX3,MX1,MX2>();
          break;
  #endif
      }
      break;
    case ROE:
      switch(dir) {
        case IDIR:
          this->RoeHD<IDIR,MX1,MX2,MX3>();
          break;
  #if DIMENSIONS >= 2
        case JDIR:
          this->RoeHD<JDIR,MX2,MX1,MX3>();
          break;
  #endif
  #if DIMENSIONS == 3
        case KDIR:
          this->RoeHD<KDIR,MX3,MX1,MX2>();
          break;
  #endif
      }
      break;
#endif // MHD
    default: // do nothing
      break;
  }

  idfx::popRegion();
}
