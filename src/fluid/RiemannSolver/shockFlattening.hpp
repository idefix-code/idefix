// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/fluid/RiemannSolver/shockFlattening.hpp
//
// Last modified : 10/2026
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2023)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_RIEMANNSOLVER_SHOCKFLATTENING_HPP_
#define FLUID_RIEMANNSOLVER_SHOCKFLATTENING_HPP_

#include <functional>

#include "fluid.hpp"
#include "../physics.hpp"
#include "eos.hpp"
#include "shockFlag.hpp"

// Forward declaration of the functor performing the actual shock detection
template<typename Phys>
struct ShockFlattening_FindShockFunctor;

template<typename Phys>
class ShockFlattening {
  // The find-shock functor reads the private members of this class
  template<typename P>
  friend struct ShockFlattening_FindShockFunctor;

 public:
  ShockFlattening(Fluid<Phys>*, real);
  ShockFlattening() {}

  // Compute the shock flag array from the current primitive state
  void FindShock();

  // Enroll a user-defined function flagging additional cells for shock flattening.
  // The function is called on the host at the end of FindShock() and receives the
  // datablock, the current time and the flag array (see UserShockFunc).
  template<typename Func>
  void EnrollUserShockFlag(Func myFunc);

  bool IsActive() const {return isActive;}
  real Smoothing() const {return smoothing;}

  // Shock flag array. Kept public so that it can be passed to (and modified by)
  // user-enrolled shock flag functions.
  IdefixArray3D<FlagShock> flagArray;

 private:
  Fluid<Phys> *hydro{nullptr};
  bool isActive{false};
  real smoothing{0};

  bool haveUserShockFlag{false};
  UserShockFunc userShockFunc;
};

template<typename Phys>
struct ShockFlattening_FindShockFunctor {
  //*****************************************************************
  // Functor constructor
  //*****************************************************************
  explicit ShockFlattening_FindShockFunctor(ShockFlattening<Phys> *sf) {
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
      eos = *(sf->hydro->eos.get());
    }
  }

  //*****************************************************************
  // Functor variables
  //*****************************************************************
  real smoothing;
  IdefixArray3D<FlagShock> flags;
  IdefixArray4D<real> Vc;
  #if GEOMETRY == CARTESIAN
    IdefixArray1D<real> dx1, dx2, dx3;
  #else
    IdefixArray3D<real> Ax1, Ax2, Ax3, dV;
  #endif

  EquationOfStateDevice eos;

  //*****************************************************************
  // Divergence of the velocity field at cell (k,j,i)
  //*****************************************************************
  KOKKOS_INLINE_FUNCTION real ComputeDivV(const int k, const int j, const int i) const {
    #if GEOMETRY == CARTESIAN
      return D_EXPAND(   (Vc(VX1,k,j,i+1) - Vc(VX1,k,j,i-1))/(2*dx1(i)) ,
                       + (Vc(VX2,k,j+1,i) - Vc(VX2,k,j-1,i))/(2*dx2(j)) ,
                       + (Vc(VX3,k+1,j,i) - Vc(VX3,k-1,j,i))/(2*dx3(k)) );
    #else
      real divV = D_EXPAND(   (Vc(VX1,k,j,i+1) + Vc(VX1,k,j,i))*Ax1(k,j,i+1) -
                              (Vc(VX1,k,j,i-1) + Vc(VX1,k,j,i))*Ax1(k,j,i)     ,
                            + (Vc(VX2,k,j+1,i) + Vc(VX2,k,j,i))*Ax2(k,j+1,i) -
                              (Vc(VX2,k,j-1,i) + Vc(VX2,k,j,i))*Ax2(k,j,i)     ,
                            + (Vc(VX3,k+1,j,i) + Vc(VX3,k,j,i))*Ax3(k+1,j,i) -
                              (Vc(VX3,k-1,j,i) + Vc(VX3,k,j,i))*Ax3(k,j,i)     );
      return 0.5*divV/dV(k,j,i);
    #endif
  }

  //*****************************************************************
  // Minimum pressure and pressure gradient magnitude around cell (k,j,i).
  // Only meaningful when Phys::pressure or Phys::isothermal.
  //*****************************************************************
  KOKKOS_INLINE_FUNCTION void ComputePressureCriterion(const int k, const int j, const int i,
                                                       real &pmin, real &gradP) const {
    if constexpr(Phys::isothermal) {
      real cs = eos.GetWaveSpeed(k,j,i);
      pmin = Vc(RHO,k,j,i)*cs*cs;
      cs = eos.GetWaveSpeed(k,j,i+1);
      real pR = Vc(RHO,k,j,i+1)*cs*cs;
      cs = eos.GetWaveSpeed(k,j,i-1);
      real pL = Vc(RHO,k,j,i-1)*cs*cs;

      pmin = FMIN(pmin,pL);
      pmin = FMIN(pmin,pR);

      gradP = FABS(pR - pL);
      #if DIMENSIONS >= 2
        cs = eos.GetWaveSpeed(k,j+1,i);
        pR = Vc(RHO,k,j+1,i)*cs*cs;
        cs = eos.GetWaveSpeed(k,j-1,i);
        pL = Vc(RHO,k,j-1,i)*cs*cs;

        pmin = FMIN(pmin,pL);
        pmin = FMIN(pmin,pR);

        gradP += FABS(pR - pL);
      #endif
      #if DIMENSIONS == 3
        cs = eos.GetWaveSpeed(k+1,j,i);
        pR = Vc(RHO,k+1,j,i)*cs*cs;
        cs = eos.GetWaveSpeed(k-1,j,i);
        pL = Vc(RHO,k-1,j,i)*cs*cs;

        pmin = FMIN(pmin,pL);
        pmin = FMIN(pmin,pR);

        gradP += FABS(pR - pL);
      #endif
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
  }

  //*****************************************************************
  // Functor operator
  //*****************************************************************
  KOKKOS_INLINE_FUNCTION void operator() (const int k, const int j, const int i) const {
    flags(k,j,i) = FlagShock::None;

    const real divV = ComputeDivV(k,j,i);
    if(divV < ZERO_F) {
      // Compressing flow: look for a strong pressure jump
      if constexpr(Phys::pressure || Phys::isothermal) {
        real pmin, gradP;
        ComputePressureCriterion(k,j,i,pmin,gradP);
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

  if (haveUserShockFlag) {
    userShockFunc(*this->hydro->data, this->hydro->data->t, this->flagArray);
  }
  idfx::popRegion();
}

template<typename Phys>
template<typename Func>
void ShockFlattening<Phys>::EnrollUserShockFlag(Func myFunc) {
  this->haveUserShockFlag = true;
  this->userShockFunc = myFunc;
}

#endif // FLUID_RIEMANNSOLVER_SHOCKFLATTENING_HPP_
