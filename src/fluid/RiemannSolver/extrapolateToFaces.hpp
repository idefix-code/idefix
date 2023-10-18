// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************
#ifndef FLUID_RIEMANNSOLVER_EXTRAPOLATETOFACES_HPP_
#define FLUID_RIEMANNSOLVER_EXTRAPOLATETOFACES_HPP_

#include "fluid.hpp"
#include "dataBlock.hpp"
#include "shockFlattening.hpp"
#include "slopeLimiter.hpp"

// Build a left and right extrapolation of the primitive variables along direction dir

// These functions extrapolate the cell prim vars to the faces. Definitions are as followed
//
// |       cell i-1               interface i          cell i
// |-----------------------------------|------------------------------------||
//          Vc(i-1)           PrimL(i)  PrimR(i)       Vc(i)
template<typename Phys,
         const int dir,
         const PLMLimiter limiter = PLMLimiter::VanLeer,
         const int order = ORDER>
class ExtrapolateToFaces {
  using SL = SlopeLimiter<limiter>;

 public:
  explicit ExtrapolateToFaces(RiemannSolver<Phys> *rSolver):
          Vc(rSolver->hydro->Vc),
          dx(rSolver->hydro->data->dx[dir]),
          shockFlattening(rSolver->haveShockFlattening),
          isRegularGrid(rSolver->hydro->data->mygrid->isRegularCartesian) {
            if(shockFlattening) {
              flags = rSolver->shockFlattening->flagArray;
            }
            if(!isRegularGrid) {
              ComputePLMweights(rSolver->hydro->data);
            }
  }

  void ComputePLMweights(DataBlock *data) {
    // Allocate weight arrays
    cpArray = IdefixArray1D<real>("ExtrapolateToFaces_cp",data->np_tot[dir]);
    cmArray = IdefixArray1D<real>("ExtrapolateToFaces_cm",data->np_tot[dir]);
    dpArray = IdefixArray1D<real>("ExtrapolateToFaces_dp",data->np_tot[dir]);
    dmArray = IdefixArray1D<real>("ExtrapolateToFaces_dm",data->np_tot[dir]);
    wpArray = IdefixArray1D<real>("ExtrapolateToFaces_wp",data->np_tot[dir]);
    wmArray = IdefixArray1D<real>("ExtrapolateToFaces_wm",data->np_tot[dir]);

    auto dx = data->dx[dir];
    auto xgc = data->xgc[dir];
    auto xr = data->xr[dir];
    auto wp = wpArray;
    auto wm = wmArray;
    auto cp = cpArray;
    auto cm = cmArray;
    auto dp = dpArray;
    auto dm = dmArray;

    idefix_for("ComputePLMweights",1,data->np_tot[dir]-1,
                KOKKOS_LAMBDA(const int i) {
                  wp(i) = dx(i) / (xgc(i+1) - xgc(i));
                  wm(i) = dx(i) / (xgc(i) - xgc(i-1));
                  cp(i) = (xgc(i+1) - xgc(i)) / (xr(i) - xgc(i));
                  cm(i) = (xgc(i) - xgc(i-1)) / (xgc(i) - xr(i-1));
                  dp(i) = (xr(i) - xgc(i)) / dx(i);
                  dm(i) = (xgc(i) - xr(i-1)) / dx(i);
                });
  }



  KOKKOS_FORCEINLINE_FUNCTION void ExtrapolatePrimVar(const int i,
                                                    const int j,
                                                    const int k,
                                                    real vL[], real vR[]) const {
    // 1-- Store the primitive variables on the left, right, and averaged states
    constexpr int ioffset = (dir==IDIR ? 1 : 0);
    constexpr int joffset = (dir==JDIR ? 1 : 0);
    constexpr int koffset = (dir==KDIR ? 1 : 0);

    for(int nv = 0 ; nv < Phys::nvar ; nv++) {
      if constexpr(order == 1) {
        vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset);
        vR[nv] = Vc(nv,k,j,i);
      } else if constexpr(order == 2) {
        if(isRegularGrid) {
          /////////////////////////////////////
          // Regular Grid, PLM reconstruction
          /////////////////////////////////////
          real dvm = Vc(nv,k-koffset,j-joffset,i-ioffset)
                    -Vc(nv,k-2*koffset,j-2*joffset,i-2*ioffset);
          real dvp = Vc(nv,k,j,i)-Vc(nv,k-koffset,j-joffset,i-ioffset);

          real dv;
          if(shockFlattening) {
            if(flags(k-koffset,j-joffset,i-ioffset) == FlagShock::Shock) {
              // Force slope limiter to minmod
              dv = SL::MinModLim(dvp,dvm);
            } else {
              dv = SL::PLMLim(dvp,dvm);
            }
          } else { // No shock flattening
            dv = SL::PLMLim(dvp,dvm);
          }

          vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset) + HALF_F*dv;

          dvm = dvp;
          dvp = Vc(nv,k+koffset,j+joffset,i+ioffset) - Vc(nv,k,j,i);

          if(shockFlattening) {
            if(flags(k,j,i) == FlagShock::Shock) {
              dv = SL::MinModLim(dvp,dvm);
            } else {
              dv = SL::PLMLim(dvp,dvm);
            }
          } else { // No shock flattening
            dv = SL::PLMLim(dvp,dvm);
          }

          vR[nv] = Vc(nv,k,j,i) - HALF_F*dv;
        } else {
          /////////////////////////////////////
          // Irregular Grid, PLM reconstruction
          /////////////////////////////////////
          const int index = ioffset*i + joffset*j + koffset*k;

          real dvm = Vc(nv,k-koffset,j-joffset,i-ioffset)
                    -Vc(nv,k-2*koffset,j-2*joffset,i-2*ioffset);
          real dvp = Vc(nv,k,j,i)-Vc(nv,k-koffset,j-joffset,i-ioffset);

          dvm *= wmArray(index-1);
          dvp *= wpArray(index-1);
          real cp = cpArray(index-1);
          real cm = cmArray(index-1);

          real dv;
          if(shockFlattening) {
            if(flags(k-koffset,j-joffset,i-ioffset) == FlagShock::Shock) {
              // Force slope limiter to minmod
              dv = SL::MinModLim(dvp,dvm);
            } else {
              dv = SL::PLMLim(dvp,dvm,cp,cm);
            }
          } else { // No shock flattening
            dv = SL::PLMLim(dvp,dvm,cp,cm);
          }

          vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset) + dpArray(index-1)*dv;

          dvm = Vc(nv,k,j,i)-Vc(nv,k-koffset,j-joffset,i-ioffset);
          dvp = Vc(nv,k+koffset,j+joffset,i+ioffset) - Vc(nv,k,j,i);
          dvm *= wmArray(index);
          dvp *= wpArray(index);
          cp = cpArray(index);
          cm = cmArray(index);

          if(shockFlattening) {
            if(flags(k,j,i) == FlagShock::Shock) {
              dv = SL::MinModLim(dvp,dvm);
            } else {
              dv = SL::PLMLim(dvp,dvm,cp,cm);
            }
          } else { // No shock flattening
            dv = SL::PLMLim(dvp,dvm,cp,cm);
          }
          vR[nv] = Vc(nv,k,j,i) - dmArray(index)*dv;
        } // Regular grid

      } else if constexpr(order == 3) {
          // 1D index along the chosen direction
          const int index = ioffset*i + joffset*j + koffset*k;
          real dvm = Vc(nv,k-koffset,j-joffset,i-ioffset)
                    -Vc(nv,k-2*koffset,j-2*joffset,i-2*ioffset);
          real dvp = Vc(nv,k,j,i)-Vc(nv,k-koffset,j-joffset,i-ioffset);

          // Limo3 limiter
          real dv;
          if(shockFlattening) {
            if(flags(k-koffset,j-joffset,i-ioffset) == FlagShock::Shock) {
              // Force slope limiter to minmod
              dv = SL::MinModLim(dvp,dvm);
            } else {
              dv = dvp * SL::LimO3Lim(dvp, dvm, dx(index-1));
            }
          } else { // No shock flattening
              dv = dvp * SL::LimO3Lim(dvp, dvm, dx(index-1));
          }

          vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset) + HALF_F*dv;

          // Check positivity
          if(nv==RHO) {
            // If face element is negative, revert to minmod
            if(vL[nv] <= 0.0) {
              dv = SL::MinModLim(dvp,dvm);
              vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset) + HALF_F*dv;
            }
          }
          if constexpr(Phys::pressure) {
            if(nv==PRS) {
              // If face element is negative, revert to minmod
              if(vL[nv] <= 0.0) {
                dv = SL::MinModLim(dvp,dvm);
                vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset) + HALF_F*dv;
              }
            }
          }

          dvm = dvp;
          dvp = Vc(nv,k+koffset,j+joffset,i+ioffset) - Vc(nv,k,j,i);

          // Limo3 limiter
          if(shockFlattening) {
            if(flags(k,j,i) == FlagShock::Shock) {
              // Force slope limiter to minmod
              dv = SL::MinModLim(dvp,dvm);
            } else {
              dv = dvm * SL::LimO3Lim(dvm, dvp, dx(index));
            }
          } else { // No shock flattening
            dv = dvm * SL::LimO3Lim(dvm, dvp, dx(index));
          }

          vR[nv] = Vc(nv,k,j,i) - HALF_F*dv;

          // Check positivity
          if(nv==RHO) {
            // If face element is negative, revert to vanleer
            if(vR[nv] <= 0.0) {
              dv = SL::MinModLim(dvp,dvm);
              vR[nv] = Vc(nv,k,j,i) - HALF_F*dv;
            }
          }
          if constexpr(Phys::pressure) {
            if(nv==PRS) {
              // If face element is negative, revert to vanleer
              if(vR[nv] <= 0.0) {
                dv = SL::MinModLim(dvp,dvm);
                vR[nv] = Vc(nv,k,j,i) - HALF_F*dv;
              }
            }
          }
      } else if constexpr(order == 4) {
          // Reconstruction in cell i-1
          real vm2 = Vc(nv,k-3*koffset,j-3*joffset,i-3*ioffset);;
          real vm1 = Vc(nv,k-2*koffset,j-2*joffset,i-2*ioffset);
          real v0 = Vc(nv,k-koffset,j-joffset,i-ioffset);
          real vp1 = Vc(nv,k,j,i);
          real vp2 = Vc(nv,k+koffset,j+joffset,i+ioffset);

          real vr,vl;
          SL::getPPMStates(vm2, vm1, v0, vp1, vp2, vl, vr);
          // vL= left side of current interface (i-1/2)= right side of cell i-1

          // Check positivity
          if(nv==RHO) {
            // If face element is negative, revert to vanleer
            if(vr <= 0.0) {
              real dv = SL::PLMLim(vp1-v0,v0-vm1);
              vr = v0+HALF_F*dv;
            }
          }
          if constexpr(Phys::pressure) {
            if(nv==PRS) {
              // If face element is negative, revert to vanleer
              if(vr <= 0.0) {
                real dv = SL::PLMLim(vp1-v0,v0-vm1);
                vr = v0+HALF_F*dv;
              }
            }
          }

          vL[nv] = vr;
          // Reconstruction in cell i

          vm2 = vm1;
          vm1 = v0;
          v0 = vp1;
          vp1 = vp2;
          vp2 = Vc(nv,k+2*koffset,j+2*joffset,i+2*ioffset);

          SL::getPPMStates(vm2, vm1, v0, vp1, vp2, vl, vr);

          // Check positivity
          if(nv==RHO) {
            // If face element is negative, revert to vanleer
            if(vl <= 0.0) {
              real dv = SL::PLMLim(vp1-v0,v0-vm1);
              vl = v0-HALF_F*dv;
            }
          }
          if constexpr(Phys::pressure) {
            if(nv==PRS) {
              // If face element is negative, revert to vanleer
              if(vl <= 0.0) {
                real dv = SL::PLMLim(vp1-v0,v0-vm1);
                vl = v0-HALF_F*dv;
              }
            }
          }

          vR[nv] = vl;
      }
    }
  }

  IdefixArray4D<real> Vc;
  IdefixArray1D<real> dx;
  IdefixArray3D<FlagShock> flags;

  IdefixArray1D<real> cpArray;
  IdefixArray1D<real> cmArray;
  IdefixArray1D<real> dpArray;
  IdefixArray1D<real> dmArray;
  IdefixArray1D<real> wpArray;
  IdefixArray1D<real> wmArray;

  bool isRegularGrid{true};
  bool shockFlattening{false};
};


#endif // FLUID_RIEMANNSOLVER_EXTRAPOLATETOFACES_HPP_
