// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "dataBlock.hpp"
#include "hydro.hpp"
#include "shockFlattening.hpp"

ShockFlattening::ShockFlattening(Hydro *h, real smoothing) {
  hydro = h;
  flagArray = IdefixArray3D<FlagShock>("flagArray",h->data->np_tot[KDIR],
                                              h->data->np_tot[JDIR],
                                              h->data->np_tot[IDIR]);
  this->isActive = true;
  this->smoothing = smoothing;
}

void ShockFlattening::FindShock() {
  idfx::pushRegion("ShockFlattening::FindShock");
  auto flags = flagArray;
  auto Vc = hydro->Vc;
  real smoothing = this->smoothing;

  #if GEOMETRY == CARTESIAN
    auto dx1 = hydro->data->dx[IDIR];
    auto dx2 = hydro->data->dx[JDIR];
    auto dx3 = hydro->data->dx[KDIR];
  #else
    auto Ax1 = hydro->data->A[IDIR];
    auto Ax2 = hydro->data->A[JDIR];
    auto Ax3 = hydro->data->A[KDIR];
    auto dV  = hydro->data->dV;
  #endif
  #ifdef ISOTHERMAL
    auto cs = hydro->isoSoundSpeedArray;
    HydroModuleStatus haveIsoCs = hydro->haveIsoSoundSpeed;
    real csIso = hydro->isoSoundSpeed;
  #endif

  auto beg = hydro->data->beg;
  auto end = hydro->data->end;
  idefix_for("findshocks",
             beg[KDIR], end[KDIR],
             beg[JDIR], end[JDIR],
             beg[IDIR], end[IDIR],
             KOKKOS_LAMBDA(int k, int j, int i) {
              flags(k,j,i) = FlagShock::None;
              #if GEOMETRY == CARTESIAN
                real divV = D_EXPAND(   (Vc(VX1,k,j,i+1) - Vc(VX1,k,j,i-1))/(2*dx1(i)) ,
                                      + (Vc(VX2,k,j+1,i) - Vc(VX2,k,j-1,i))/(2*dx2(j)) ,
                                      + (Vc(VX3,k+1,j,i) - Vc(VX3,k-1,j,i))/(2*dx3(k)) );
              #else
                real divV = D_EXPAND(   (Vc(VX1,k,j,i+1) + Vc(VX1,k,j,i))*Ax1(k,j,i+1) -
                                        (Vc(VX1,k,j,i-1) + Vc(VX1,k,j,i))*Ax1(k,j,i)     ,
                                      + (Vc(VX2,k,j+1,i) + Vc(VX2,k,j,i))*Ax2(k,j+1,i) -
                                        (Vc(VX2,k,j-1,i) + Vc(VX2,k,j,i))*Ax2(k,j,i)     ,
                                      + (Vc(VX3,k+1,j,i) + Vc(VX3,k,j,i))*Ax3(k+1,j,i) -
                                        (Vc(VX3,k-1,j,i) + Vc(VX3,k,j,i))*Ax3(k,j,i)     );
                divV = 0.5*divV/dV(k,j,i);
              #endif

              if(divV<ZERO_F) {
                #ifdef ISOTHERMAL
                  real pmin, gradP;
                  if(haveIsoCs == UserDefFunction) {
                    pmin = Vc(RHO,k,j,i)*cs(k,j,i)*cs(k,j,i);
                    pmin = FMIN(pmin,Vc(RHO,k,j,i+1)*cs(k,j,i+1)*cs(k,j,i+1));
                    pmin = FMIN(pmin,Vc(RHO,k,j,i-1)*cs(k,j,i-1)*cs(k,j,i-1));
                    gradP = FABS(Vc(RHO,k,j,i+1)*cs(k,j,i+1)*cs(k,j,i+1) -
                                 Vc(RHO,k,j,i-1)*cs(k,j,i-1)*cs(k,j,i-1));
                    #if DIMENSIONS >= 2
                      pmin = FMIN(pmin,Vc(RHO,k,j+1,i)*cs(k,j+1,i)*cs(k,j+1,i));
                      pmin = FMIN(pmin,Vc(RHO,k,j-1,i)*cs(k,j-1,i)*cs(k,j-1,i));
                      gradP += FABS(Vc(RHO,k,j+1,i)*cs(k,j+1,i)*cs(k,j+1,i) -
                                    Vc(RHO,k,j-1,i)*cs(k,j-1,i)*cs(k,j-1,i));
                    #endif
                    #if DIMENSIONS == 3
                      pmin = FMIN(pmin,Vc(RHO,k+1,j,i)*cs(k+1,j,i)*cs(k+1,j,i));
                      pmin = FMIN(pmin,Vc(RHO,k-1,j,i)*cs(k-1,j,i)*cs(k-1,j,i));
                      gradP += FABS(Vc(RHO,k+1,j,i)*cs(k+1,j,i)*cs(k+1,j,i) -
                                    Vc(RHO,k-1,j,i)*cs(k-1,j,i)*cs(k-1,j,i));
                    #endif
                  } else {
                    pmin = Vc(RHO,k,j,i)*csIso*csIso;
                    pmin = FMIN(pmin,Vc(RHO,k,j,i+1)*csIso*csIso);
                    pmin = FMIN(pmin,Vc(RHO,k,j,i-1)*csIso*csIso);
                    gradP = FABS(Vc(RHO,k,j,i+1)*csIso*csIso -
                                 Vc(RHO,k,j,i-1)*csIso*csIso);
                    #if DIMENSIONS >= 2
                      pmin = FMIN(pmin,Vc(RHO,k,j+1,i)*csIso*csIso);
                      pmin = FMIN(pmin,Vc(RHO,k,j-1,i)*csIso*csIso);
                      gradP += FABS(Vc(RHO,k,j+1,i)*csIso*csIso -
                                    Vc(RHO,k,j-1,i)*csIso*csIso);
                    #endif
                    #if DIMENSIONS == 3
                      pmin = FMIN(pmin,Vc(RHO,k+1,j,i)*csIso*csIso);
                      pmin = FMIN(pmin,Vc(RHO,k-1,j,i)*csIso*csIso);
                      gradP += FABS(Vc(RHO,k+1,j,i)*csIso*csIso -
                                    Vc(RHO,k-1,j,i)*csIso*csIso);
                    #endif
                  }
                #else // ISOTHERMAL
                  real pmin = Vc(PRS,k,j,i);
                  pmin = FMIN(pmin,Vc(PRS,k,j,i+1));
                  pmin = FMIN(pmin,Vc(PRS,k,j,i-1));
                  real gradP = FABS(Vc(PRS,k,j,i+1) - Vc(PRS,k,j,i-1));
                  #if DIMENSIONS >= 2
                    pmin = FMIN(pmin,Vc(PRS,k,j+1,i));
                    pmin = FMIN(pmin,Vc(PRS,k,j-1,i));
                    gradP += FABS(Vc(PRS,k,j+1,i) - Vc(PRS,k,j-1,i));
                  #endif
                  #if DIMENSIONS == 3
                    pmin = FMIN(pmin,Vc(PRS,k+1,j,i));
                    pmin = FMIN(pmin,Vc(PRS,k-1,j,i));
                    gradP += FABS(Vc(PRS,k+1,j,i) - Vc(PRS,k-1,j,i));
                  #endif
                #endif // ISOTHERMAL
                if(gradP > smoothing*pmin) {
                  flags(k,j,i) = FlagShock::Shock;
                }
              }
            });
  idfx::popRegion();
}
