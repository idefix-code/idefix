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
void Hydro::CalcRiemannFlux(DataBlock & data, int dir, const real t) {
  idfx::pushRegion("Hydro::CalcRiemannFlux");

  if(haveHall == UserDefFunction && dir == IDIR) {
    if(hallDiffusivityFunc)
      hallDiffusivityFunc(data, t, data.xHall);
    else
      IDEFIX_ERROR("No user-defined Hall diffusivity function has been enrolled");
  }

  switch (mySolver) {
#if MHD == YES
    case TVDLF:
      switch(dir) {
        case IDIR:
          TvdlfMHD<IDIR,ARG_EXPAND(MX1,MX2,MX3),
                   ARG_EXPAND(BX1,BX2,BX3)>(data, this->gamma, this->C2Iso);
          break;
#if DIMENSIONS >= 2
        case JDIR:
          TvdlfMHD<JDIR,ARG_EXPAND(MX2,MX1,MX3),
                   ARG_EXPAND(BX2,BX1,BX3)>(data, this->gamma, this->C2Iso);
          break;
#endif
#if DIMENSIONS == 3
        case KDIR:
          TvdlfMHD<KDIR,ARG_EXPAND(MX3,MX1,MX2),
                   ARG_EXPAND(BX3,BX1,BX2)>(data, this->gamma, this->C2Iso);
          break;
#endif
      }
      break;
    case HLL:
      switch(dir) {
        case IDIR:
          HllMHD<IDIR,ARG_EXPAND(MX1,MX2,MX3),
                 ARG_EXPAND(BX1,BX2,BX3)>(data, this->gamma, this->C2Iso,
                                               this->haveHall, this->xH);
          break;
#if DIMENSIONS >= 2
        case JDIR:
          HllMHD<JDIR,ARG_EXPAND(MX2,MX1,MX3),
                 ARG_EXPAND(BX2,BX1,BX3)>(data, this->gamma, this->C2Iso,
                                               this->haveHall, this->xH);
          break;
#endif
#if DIMENSIONS == 3
        case KDIR:
          HllMHD<KDIR,ARG_EXPAND(MX3,MX1,MX2),
                 ARG_EXPAND(BX3,BX1,BX2)>(data, this->gamma, this->C2Iso,
                                               this->haveHall, this->xH);
          break;
#endif
      }
      break;
    case HLLD:
      switch(dir) {
        case IDIR:
          HlldMHD<IDIR,ARG_EXPAND(MX1,MX2,MX3),
                  ARG_EXPAND(BX1,BX2,BX3)>(data, this->gamma, this->C2Iso);
          break;
#if DIMENSIONS >= 2
        case JDIR:
          HlldMHD<JDIR,ARG_EXPAND(MX2,MX1,MX3),
                  ARG_EXPAND(BX2,BX1,BX3)>(data, this->gamma, this->C2Iso);
          break;
#endif
#if DIMENSIONS == 3
        case KDIR:
          HlldMHD<KDIR,ARG_EXPAND(MX3,MX1,MX2),
                  ARG_EXPAND(BX3,BX1,BX2)>(data, this->gamma, this->C2Iso);
          break;
#endif
      }
      break;
    case ROE:
      switch(dir) {
        case IDIR:
          RoeMHD<IDIR,ARG_EXPAND(MX1,MX2,MX3),
                 ARG_EXPAND(BX1,BX2,BX3)>(data, this->gamma, this->C2Iso);
          break;
#if DIMENSIONS >= 2
        case JDIR:
          RoeMHD<JDIR,ARG_EXPAND(MX2,MX1,MX3),
                 ARG_EXPAND(BX2,BX1,BX3)>(data, this->gamma, this->C2Iso);
          break;
#endif
#if DIMENSIONS == 3
        case KDIR:
          RoeMHD<KDIR,ARG_EXPAND(MX3,MX1,MX2),
                 ARG_EXPAND(BX3,BX1,BX2)>(data, this->gamma, this->C2Iso);
          break;
#endif
      }
      break;
#else
    case TVDLF: 
      switch(dir) {
        case IDIR:
          TvdlfHD<IDIR,MX1,MX2,MX3>(data, this->gamma, this->C2Iso);
          break;
#if DIMENSIONS >= 2
        case JDIR:
          TvdlfHD<JDIR,MX2,MX1,MX3>(data, this->gamma, this->C2Iso);
          break;
#endif
#if DIMENSIONS == 3
        case KDIR:
          TvdlfHD<KDIR,MX3,MX1,MX2>(data, this->gamma, this->C2Iso);
          break;
#endif
      }
      break;
    case HLL:
      switch(dir) {
        case IDIR:
          HllHD<IDIR,MX1,MX2,MX3>(data, this->gamma, this->C2Iso);
          break;
#if DIMENSIONS >= 2
        case JDIR:
          HllHD<JDIR,MX2,MX1,MX3>(data, this->gamma, this->C2Iso);
          break;
#endif
#if DIMENSIONS == 3
        case KDIR:
          HllHD<KDIR,MX3,MX1,MX2>(data, this->gamma, this->C2Iso);
          break;
#endif
      }
      break;
    case HLLC:
      switch(dir) {
        case IDIR:
          HllcHD<IDIR,MX1,MX2,MX3>(data, this->gamma, this->C2Iso);
          break;
#if DIMENSIONS >= 2
        case JDIR:
          HllcHD<JDIR,MX2,MX1,MX3>(data, this->gamma, this->C2Iso);
          break;
#endif
#if DIMENSIONS == 3
        case KDIR:
          HllcHD<KDIR,MX3,MX1,MX2>(data, this->gamma, this->C2Iso);
          break;
#endif
      }
      break;
    case ROE:
      switch(dir) {
        case IDIR:
          RoeHD<IDIR,MX1,MX2,MX3>(data, this->gamma, this->C2Iso);
          break;
#if DIMENSIONS >= 2
        case JDIR:
          RoeHD<JDIR,MX2,MX1,MX3>(data, this->gamma, this->C2Iso);
          break;
#endif
#if DIMENSIONS == 3
        case KDIR:
          RoeHD<KDIR,MX3,MX1,MX2>(data, this->gamma, this->C2Iso);
          break;
#endif
      }
      break;
#endif
    default: // do nothing
      break;
  }

  idfx::popRegion();
}
