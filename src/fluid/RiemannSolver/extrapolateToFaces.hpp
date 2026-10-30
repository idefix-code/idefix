// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/fluid/RiemannSolver/extrapolateToFaces.hpp
//
// Last modified : 10/2026
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2023)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************
#ifndef FLUID_RIEMANNSOLVER_EXTRAPOLATETOFACES_HPP_
#define FLUID_RIEMANNSOLVER_EXTRAPOLATETOFACES_HPP_

#include "dataBlock.hpp"
#include "shockFlag.hpp"
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

  // Index offsets along the reconstruction direction
  static constexpr int ioffset = (dir == IDIR ? 1 : 0);
  static constexpr int joffset = (dir == JDIR ? 1 : 0);
  static constexpr int koffset = (dir == KDIR ? 1 : 0);

 public:
  // data: datablock holding the grid information
  // VcArray: primitive variables in cell centers
  // shockFlags: shock flag array (only used when hasShockFlattening is true)
  explicit ExtrapolateToFaces(DataBlock *data,
                              IdefixArray4D<real> VcArray,
                              IdefixArray3D<FlagShock> shockFlags = {},
                              bool hasShockFlattening = false):
          Vc(VcArray),
          dx(data->dx[dir]),
          flags(shockFlags),
          isRegularGrid(data->mygrid->isRegularCartesian),
          shockFlattening(hasShockFlattening) {
    if(!isRegularGrid) {
      ComputeGridWeights(data);
    }
  }

  // Extrapolate the primitive variables to the left (vL) and right (vR) states
  // of the interface located between cell i-1 and cell i
  KOKKOS_FORCEINLINE_FUNCTION void ExtrapolatePrimVar(const int i,
                                                      const int j,
                                                      const int k,
                                                      real vL[], real vR[]) const {
    if constexpr(order == 1) {
      ExtrapolateOrder1(i,j,k,vL,vR);
    } else if constexpr(order == 2) {
      if(isRegularGrid) {
        ExtrapolateOrder2Regular(i,j,k,vL,vR);
      } else {
        ExtrapolateOrder2Irregular(i,j,k,vL,vR);
      }
    } else if constexpr(order == 3) {
      ExtrapolateOrder3(i,j,k,vL,vR);
    } else if constexpr(order == 4) {
      ExtrapolateOrder4(i,j,k,vL,vR);
    }
  }

 private:
  // PLM weights for irregular grids (trivially copyable, hence device-safe)
  struct GridWeights {
    IdefixArray1D<real> cp;   ///< (xgc(i+1)-xgc(i)) / (xr(i)-xgc(i))
    IdefixArray1D<real> cm;   ///< (xgc(i)-xgc(i-1)) / (xgc(i)-xr(i-1))
    IdefixArray1D<real> dp;   ///< (xr(i)-xgc(i)) / dx(i)
    IdefixArray1D<real> dm;   ///< (xgc(i)-xr(i-1)) / dx(i)
    IdefixArray1D<real> wp;   ///< dx(i) / (xgc(i+1)-xgc(i))
    IdefixArray1D<real> wm;   ///< dx(i) / (xgc(i)-xgc(i-1))
  };

  // Allocate and compute the PLM weights along direction dir (irregular grids only)
  void ComputeGridWeights(DataBlock *data) {
    gridWeights.cp = IdefixArray1D<real>("ExtrapolateToFaces_cp",data->np_tot[dir]);
    gridWeights.cm = IdefixArray1D<real>("ExtrapolateToFaces_cm",data->np_tot[dir]);
    gridWeights.dp = IdefixArray1D<real>("ExtrapolateToFaces_dp",data->np_tot[dir]);
    gridWeights.dm = IdefixArray1D<real>("ExtrapolateToFaces_dm",data->np_tot[dir]);
    gridWeights.wp = IdefixArray1D<real>("ExtrapolateToFaces_wp",data->np_tot[dir]);
    gridWeights.wm = IdefixArray1D<real>("ExtrapolateToFaces_wm",data->np_tot[dir]);

    auto dxGrid = data->dx[dir];
    auto xgc = data->xgc[dir];
    auto xr = data->xr[dir];
    auto wp = gridWeights.wp;
    auto wm = gridWeights.wm;
    auto cp = gridWeights.cp;
    auto cm = gridWeights.cm;
    auto dp = gridWeights.dp;
    auto dm = gridWeights.dm;

    idefix_for("ComputeGridWeights",1,data->np_tot[dir]-1,
               KOKKOS_LAMBDA(const int i) {
                 wp(i) = dxGrid(i) / (xgc(i+1) - xgc(i));
                 wm(i) = dxGrid(i) / (xgc(i) - xgc(i-1));
                 cp(i) = (xgc(i+1) - xgc(i)) / (xr(i) - xgc(i));
                 cm(i) = (xgc(i) - xgc(i-1)) / (xgc(i) - xr(i-1));
                 dp(i) = (xr(i) - xgc(i)) / dxGrid(i);
                 dm(i) = (xgc(i) - xr(i-1)) / dxGrid(i);
               });
  }

  // Shock flag at cell (k,j,i), or FlagShock::None when shock flattening is disabled
  KOKKOS_FORCEINLINE_FUNCTION FlagShock GetFlag(const int k, const int j, const int i) const {
    return shockFlattening ? flags(k,j,i) : FlagShock::None;
  }

  // PLM slope with a minmod fallback in shocked cells (regular grid)
  KOKKOS_FORCEINLINE_FUNCTION real LimitSlope(const real dvp, const real dvm,
                                               const FlagShock flag) const {
    if(shockFlattening && flag == FlagShock::Shock) {
      return SL::MinModLim(dvp,dvm);
    }
    return SL::PLMLim(dvp,dvm);
  }

  // PLM slope with a minmod fallback in shocked cells (irregular grid)
  KOKKOS_FORCEINLINE_FUNCTION real LimitSlope(const real dvp, const real dvm,
                                               const FlagShock flag,
                                               const real cp, const real cm) const {
    if(shockFlattening && flag == FlagShock::Shock) {
      return SL::MinModLim(dvp,dvm);
    }
    return SL::PLMLim(dvp,dvm,cp,cm);
  }

  // LimO3 slope with a minmod fallback in shocked cells.
  // dvU multiplies the limiter, (dvA,dvB) are the limiter arguments.
  KOKKOS_FORCEINLINE_FUNCTION real LimitSlopeO3(const real dvU, const real dvA, const real dvB,
                                                 const FlagShock flag, const real dxm) const {
    if(shockFlattening && flag == FlagShock::Shock) {
      return SL::MinModLim(dvA,dvB);
    }
    return dvU * SL::LimO3Lim(dvA,dvB,dxm);
  }

  // Enforce positivity of the density and pressure face values by reverting to the
  // fallback slope dvFallback (sign = +1 for a left state, -1 for a right state)
  KOKKOS_FORCEINLINE_FUNCTION void EnforcePositivity(const int nv, real &vFace,
                                                      const real v0, const real dvFallback,
                                                      const real sign) const {
    bool enforce = (nv == RHO);
    if constexpr(Phys::pressure) {
      enforce = enforce || (nv == PRS);
    }
    if(enforce && vFace <= 0.0) {
      vFace = v0 + sign*HALF_F*dvFallback;
    }
  }

  // First order: piecewise constant reconstruction
  KOKKOS_FORCEINLINE_FUNCTION void ExtrapolateOrder1(const int i, const int j, const int k,
                                                      real vL[], real vR[]) const {
    for(int nv = 0 ; nv < Phys::nvar ; nv++) {
      vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset);
      vR[nv] = Vc(nv,k,j,i);
    }
  }

  // Second order (PLM) reconstruction on a regular grid
  KOKKOS_FORCEINLINE_FUNCTION void ExtrapolateOrder2Regular(const int i, const int j, const int k,
                                                             real vL[], real vR[]) const {
    for(int nv = 0 ; nv < Phys::nvar ; nv++) {
      // Left state: reconstruction in cell i-1
      real dvm = Vc(nv,k-koffset,j-joffset,i-ioffset)
                -Vc(nv,k-2*koffset,j-2*joffset,i-2*ioffset);
      real dvp = Vc(nv,k,j,i)-Vc(nv,k-koffset,j-joffset,i-ioffset);
      real dv = LimitSlope(dvp,dvm,GetFlag(k-koffset,j-joffset,i-ioffset));
      vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset) + HALF_F*dv;

      // Right state: reconstruction in cell i
      dvm = dvp;
      dvp = Vc(nv,k+koffset,j+joffset,i+ioffset) - Vc(nv,k,j,i);
      dv = LimitSlope(dvp,dvm,GetFlag(k,j,i));
      vR[nv] = Vc(nv,k,j,i) - HALF_F*dv;
    }
  }

  // Second order (PLM) reconstruction on an irregular grid
  KOKKOS_FORCEINLINE_FUNCTION void ExtrapolateOrder2Irregular(const int i, const int j, const int k,
                                                               real vL[], real vR[]) const {
    // 1D index along the chosen direction
    const int index = ioffset*i + joffset*j + koffset*k;

    for(int nv = 0 ; nv < Phys::nvar ; nv++) {
      // Left state: reconstruction in cell i-1
      real dvm = Vc(nv,k-koffset,j-joffset,i-ioffset)
                -Vc(nv,k-2*koffset,j-2*joffset,i-2*ioffset);
      real dvp = Vc(nv,k,j,i)-Vc(nv,k-koffset,j-joffset,i-ioffset);
      dvm *= gridWeights.wm(index-1);
      dvp *= gridWeights.wp(index-1);
      real dv = LimitSlope(dvp,dvm,GetFlag(k-koffset,j-joffset,i-ioffset),
                           gridWeights.cp(index-1),gridWeights.cm(index-1));
      vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset) + gridWeights.dp(index-1)*dv;

      // Right state: reconstruction in cell i
      dvm = Vc(nv,k,j,i)-Vc(nv,k-koffset,j-joffset,i-ioffset);
      dvp = Vc(nv,k+koffset,j+joffset,i+ioffset) - Vc(nv,k,j,i);
      dvm *= gridWeights.wm(index);
      dvp *= gridWeights.wp(index);
      dv = LimitSlope(dvp,dvm,GetFlag(k,j,i),
                      gridWeights.cp(index),gridWeights.cm(index));
      vR[nv] = Vc(nv,k,j,i) - gridWeights.dm(index)*dv;
    }
  }

  // Third order (LimO3) reconstruction
  KOKKOS_FORCEINLINE_FUNCTION void ExtrapolateOrder3(const int i, const int j, const int k,
                                                      real vL[], real vR[]) const {
    // 1D index along the chosen direction
    const int index = ioffset*i + joffset*j + koffset*k;

    for(int nv = 0 ; nv < Phys::nvar ; nv++) {
      // Left state: reconstruction in cell i-1
      real dvm = Vc(nv,k-koffset,j-joffset,i-ioffset)
                -Vc(nv,k-2*koffset,j-2*joffset,i-2*ioffset);
      real dvp = Vc(nv,k,j,i)-Vc(nv,k-koffset,j-joffset,i-ioffset);
      real dv = LimitSlopeO3(dvp,dvp,dvm,GetFlag(k-koffset,j-joffset,i-ioffset),dx(index-1));
      vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset) + HALF_F*dv;
      // Check positivity of density and pressure
      EnforcePositivity(nv,vL[nv],Vc(nv,k-koffset,j-joffset,i-ioffset),
                        SL::MinModLim(dvp,dvm),ONE_F);

      // Right state: reconstruction in cell i
      dvm = dvp;
      dvp = Vc(nv,k+koffset,j+joffset,i+ioffset) - Vc(nv,k,j,i);
      dv = LimitSlopeO3(dvm,dvm,dvp,GetFlag(k,j,i),dx(index));
      vR[nv] = Vc(nv,k,j,i) - HALF_F*dv;
      // Check positivity of density and pressure
      EnforcePositivity(nv,vR[nv],Vc(nv,k,j,i),SL::MinModLim(dvp,dvm),-ONE_F);
    }
  }

  // Fourth order (PPM) reconstruction
  KOKKOS_FORCEINLINE_FUNCTION void ExtrapolateOrder4(const int i, const int j, const int k,
                                                      real vL[], real vR[]) const {
    for(int nv = 0 ; nv < Phys::nvar ; nv++) {
      // Reconstruction in cell i-1
      real vm2 = Vc(nv,k-3*koffset,j-3*joffset,i-3*ioffset);
      real vm1 = Vc(nv,k-2*koffset,j-2*joffset,i-2*ioffset);
      real v0  = Vc(nv,k-koffset,j-joffset,i-ioffset);
      real vp1 = Vc(nv,k,j,i);
      real vp2 = Vc(nv,k+koffset,j+joffset,i+ioffset);

      real vr, vl;
      SL::getPPMStates(vm2,vm1,v0,vp1,vp2,vl,vr);
      // vL is the left side of the current interface, i.e. the right side of cell i-1
      EnforcePositivity(nv,vr,v0,SL::PLMLim(vp1-v0,v0-vm1),ONE_F);
      vL[nv] = vr;

      // Reconstruction in cell i
      vm2 = vm1;
      vm1 = v0;
      v0 = vp1;
      vp1 = vp2;
      vp2 = Vc(nv,k+2*koffset,j+2*joffset,i+2*ioffset);

      SL::getPPMStates(vm2,vm1,v0,vp1,vp2,vl,vr);
      // vR is the right side of the current interface, i.e. the left side of cell i
      EnforcePositivity(nv,vl,v0,SL::PLMLim(vp1-v0,v0-vm1),-ONE_F);
      vR[nv] = vl;
    }
  }

  IdefixArray4D<real> Vc;            ///< primitive variables in cell centers
  IdefixArray1D<real> dx;            ///< cell width along direction dir
  IdefixArray3D<FlagShock> flags;    ///< shock flags (empty when shock flattening is disabled)
  GridWeights gridWeights;           ///< PLM weights (irregular grids only)

  bool isRegularGrid{true};
  bool shockFlattening{false};
};


#endif // FLUID_RIEMANNSOLVER_EXTRAPOLATETOFACES_HPP_
