// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_RIEMANNSOLVER_SHOCKFLATTENING_HPP_
#define FLUID_RIEMANNSOLVER_SHOCKFLATTENING_HPP_

#include "fluid.hpp"
#include "../physics.hpp"

enum class FlagShock{None, Shock};

template<typename Phys>
class ShockFlattening {
 public:
  ShockFlattening(Fluid<Phys>*, real);
  ShockFlattening() {}

  void FindShock();

  Fluid<Phys> *hydro;
  IdefixArray3D<FlagShock> flagArray;
  bool isActive{false};
  real smoothing{0};
};

template<typename Phys>
struct ShockFlattening_FindShockFunctor {
  //*****************************************************************
  // Functor constructor
  //*****************************************************************
  explicit ShockFlattening_FindShockFunctor(ShockFlattening<Phys>  *sf) {
    flags = sf->flagArray;
    Vc = sf->hydro->Vc;
    smoothing = sf->smoothing;

    #if GEOMETRY == CARTESIAN
      dx1 = sf->hydro->data->dx[IDIR];
      dx2 = sf->hydro->data->dx[JDIR];
      dx3 = sf->hydro->data->dx[KDIR];
    #else
      Ax1 = sf->hydro->data->A[IDIR];
      Ax2 = sf->hydro->data->A[JDIR];
      Ax3 = sf->hydro->data->A[KDIR];
      dV  = sf->hydro->data->dV;
    #endif
    if constexpr(Phys::isothermal) {
      cs = sf->hydro->isoSoundSpeedArray;
      haveIsoCs = sf->hydro->haveIsoSoundSpeed;
      csIso = sf->hydro->isoSoundSpeed;
    }
  }
  //*****************************************************************
  // Functor Variables
  //*****************************************************************
  real smoothing;
  IdefixArray3D<FlagShock> flags;
  IdefixArray4D<real> Vc;
  #if GEOMETRY == CARTESIAN
    IdefixArray1D<real> dx1, dx2, dx3;
  #else
    IdefixArray3D<real> Ax1,Ax2,Ax3, dV;
  #endif

  IdefixArray3D<real> cs;
  HydroModuleStatus haveIsoCs;;
  real csIso;

    //*****************************************************************
  // Functor Operator
  //*****************************************************************
  KOKKOS_INLINE_FUNCTION void operator() (const int k, const int j,  const int i) const {
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
      real pmin, gradP;
      if constexpr(Phys::isothermal) {
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
      } else if constexpr(Phys::pressure) {
        pmin = Vc(PRS,k,j,i);
        pmin = FMIN(pmin,Vc(PRS,k,j,i+1));
        pmin = FMIN(pmin,Vc(PRS,k,j,i-1));
        gradP = FABS(Vc(PRS,k,j,i+1) - Vc(PRS,k,j,i-1));
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
      }
      if constexpr(Phys::pressure || Phys::isothermal) {
        if(gradP > smoothing*pmin) {
          flags(k,j,i) = FlagShock::Shock;
        }
      }
    }
  }
};


template<typename Phys>
ShockFlattening<Phys>::ShockFlattening(Fluid<Phys> *h, real smoothing) {
  hydro = h;
  flagArray = IdefixArray3D<FlagShock>("flagArray",h->data->np_tot[KDIR],
                                              h->data->np_tot[JDIR],
                                              h->data->np_tot[IDIR]);
  this->isActive = true;
  this->smoothing = smoothing;
}

template<typename Phys>
void ShockFlattening<Phys>::FindShock() {
  idfx::pushRegion("ShockFlattening::FindShock");

  auto func = ShockFlattening_FindShockFunctor<Phys>(this);

  auto beg = hydro->data->beg;
  auto end = hydro->data->end;
  idefix_for("findshocks",
             beg[KDIR]-KOFFSET, end[KDIR]+KOFFSET,
             beg[JDIR]-JOFFSET, end[JDIR]+JOFFSET,
             beg[IDIR]-IOFFSET, end[IDIR]+IOFFSET,
             func);
  idfx::popRegion();
}


#endif // FLUID_RIEMANNSOLVER_SHOCKFLATTENING_HPP_
