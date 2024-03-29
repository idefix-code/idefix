// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************
#ifndef FLUID_CALCCURRENT_HPP_
#define FLUID_CALCCURRENT_HPP_

#include "fluid.hpp"
#include "dataBlock.hpp"

// Compute the electrical current on faces
template <typename Phys>
void Fluid<Phys>::CalcCurrent() {
  idfx::pushRegion("Fluid::CalcCurrent");
  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Vs = this->Vs;
  IdefixArray4D<real> J = this->J;

  IdefixArray1D<real> dx1 = data->dx[IDIR];
  IdefixArray1D<real> dx2 = data->dx[JDIR];
  IdefixArray1D<real> dx3 = data->dx[KDIR];

  IdefixArray1D<real> r = data->x[IDIR];
  IdefixArray1D<real> rm = data->xl[IDIR];
  IdefixArray1D<real> th = data->x[JDIR];

  #if GEOMETRY == SPHERICAL
  IdefixArray1D<real> sinx2 = data->sinx2;
  IdefixArray1D<real> sinx2m = data->sinx2m;
  #endif

  bool haveGridCoarsening = data->haveGridCoarsening != GridCoarsening::disabled;
  bool haveGridCoarseningX1 = data->coarseningDirection[IDIR];
  bool haveGridCoarseningX2 = data->coarseningDirection[JDIR];
  bool haveGridCoarseningX3 = data->coarseningDirection[KDIR];
  auto coarseningLevelX1 = data->coarseningLevel[IDIR];
  auto coarseningLevelX2 = data->coarseningLevel[JDIR];
  auto coarseningLevelX3 = data->coarseningLevel[KDIR];

  idefix_for("CalcCurrent",
             KOFFSET,data->np_tot[KDIR],
             JOFFSET,data->np_tot[JDIR],
             IOFFSET,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      [[maybe_unused]] real Bx1_000 = ZERO_F, Bx1_0m0 = ZERO_F, Bx1_00m = ZERO_F;
      [[maybe_unused]] real Bx2_000 = ZERO_F, Bx2_m00 = ZERO_F, Bx2_00m = ZERO_F;
      [[maybe_unused]] real Bx3_000 = ZERO_F, Bx3_m00 = ZERO_F, Bx3_0m0 = ZERO_F;

      [[maybe_unused]] real d12 = ZERO_F, d13 = ZERO_F,
                            d21 = ZERO_F, d23 = ZERO_F,
                            d31 = ZERO_F, d32 = ZERO_F;

      [[maybe_unused]] real Jx = ZERO_F, Jy = ZERO_F, Jz = ZERO_F;

      Bx1_000 = Vs(BX1s,k,j,i);

#if DIMENSIONS >= 2
      Bx1_0m0 = Vs(BX1s,k,j-1,i);
      Bx2_000 = Vs(BX2s,k,j,i);
      Bx2_m00 = Vs(BX2s,k,j,i-1);
#elif COMPONENTS >= 2   // DIMENSIONS == 1 here!
      Bx2_000 = Vc(BX2,k,j,i);
      Bx2_m00 = Vc(BX2,k,j,i-1);
  #if COMPONENTS == 3
      Bx3_000 = Vc(BX3,k,j,i);
      Bx3_m00 = Vc(BX3,k,j,i-1);
  #endif
#endif

#if DIMENSIONS == 3
      Bx1_00m = Vs(BX1s,k-1,j,i);
      Bx2_00m = Vs(BX2s,k-1,j,i);
      Bx3_000 = Vs(BX3s,k,j,i);
      Bx3_m00 = Vs(BX3s,k,j,i-1);
      Bx3_0m0 = Vs(BX3s,k,j-1,i);
#elif COMPONENTS == 3  && DIMENSIONS == 2 // DIMENSIONS 2 here
      Bx3_000 = Vc(BX3,k,j,i);
      Bx3_m00 = Vc(BX3,k,j,i-1);
      Bx3_0m0 = Vc(BX3,k,j-1,i);
#endif


      // define geometrical factors
      D_EXPAND( d12 = d13 = TWO_F / (dx1(i-1)+dx1(i));  ,
                d21 = d23 = TWO_F / (dx2(j-1)+dx2(j));  ,
                d31 = d32 = TWO_F / (dx3(k-1)+dx3(k));  )

      [[maybe_unused]] real a13Bx3_000 = Bx3_000;
      [[maybe_unused]] real a13Bx3_m00 = Bx3_m00;
      [[maybe_unused]] real a23Bx3_000 = Bx3_000;
      [[maybe_unused]] real a23Bx3_0m0 = Bx3_0m0;
      [[maybe_unused]] real a12Bx2_000 = Bx2_000;
      [[maybe_unused]] real a12Bx2_m00 = Bx2_m00;


#if GEOMETRY == CYLINDRICAL
      d32 = d31 = ZERO_F;
      d13 = TWO_F / (FABS(r(i))*r(i) - FABS(r(i-1))*r(i-1));
  #if COMPONENTS == 3
      a13Bx3_000 = Bx3_000 * FABS(r(i));
      a13Bx3_m00 = Bx3_m00 * FABS(r(i-1));
  #endif

#elif GEOMETRY == POLAR
      d23 /= r(i);
      d21 /= rm(i);
      d12 = TWO_F / (FABS(r(i))*r(i) - FABS(r(i-1))*r(i-1));
  #if COMPONENTS >= 2
      a12Bx2_000 = Bx2_000 * r(i);
      a12Bx2_m00 = Bx2_m00 * r(i-1);
  #endif

#elif GEOMETRY == SPHERICAL
      real s = FABS(sinx2(j));
      real sm = 0.5*(FABS(sinx2(j))+FABS(sinx2(j-1)));

      D_EXPAND(d12 /= rm(i);   d13 /= rm(i);   ,
              d21 /= rm(i);   d23 /= r(i)*sm;  ,
              d32 /= r(i)*sm; d31 /= rm(i)*s;  )

  #if COMPONENTS >= 2
      a12Bx2_000 = Bx2_000 * r(i);
      a12Bx2_m00 = Bx2_m00 * r(i-1);
  #endif
  #if COMPONENTS == 3
      a13Bx3_000 = Bx3_000 * r(i);
      a13Bx3_m00 = Bx3_m00 * r(i-1);
    #if DIMENSIONS >= 2
      a23Bx3_000 = Bx3_000 * s;
      a23Bx3_0m0 = Bx3_0m0 * FABS(sinx2(j-1));
    #endif
  #endif // COMPONENTS

#endif

    // Correct spacing for grid coarsening
    if(haveGridCoarsening) {
      if(haveGridCoarseningX1) {
        const int factor =  1 << (coarseningLevelX1(k,j) - 1);
        d12 /= factor;
        d13 /= factor;
      }

      if(haveGridCoarseningX2) {
        const int factor =  1 << (coarseningLevelX2(k,i) - 1);
        d21 /= factor;
        d23 /= factor;
      }

      if(haveGridCoarseningX3) {
        const int factor =  1 << (coarseningLevelX3(j,i) - 1);
        d31 /= factor;
        d32 /= factor;
      }
    }

    // Compute actual current

#if COMPONENTS == 3
  #if DIMENSIONS == 3
      Jx +=  - (Bx2_000 - Bx2_00m)*d32;
      Jy +=  + (Bx1_000 - Bx1_00m)*d31;
  #endif
  #if DIMENSIONS >= 2
      Jx += + (a23Bx3_000 - a23Bx3_0m0)*d23;
  #endif
      Jy += - (a13Bx3_000 - a13Bx3_m00)*d13;
#endif

#if COMPONENTS >= 2
      Jz += (a12Bx2_000 - a12Bx2_m00)*d12;
#endif

#if DIMENSIONS >= 2
      Jz += - (Bx1_000 - Bx1_0m0)*d21;
#endif

      // Store the beast
      J(IDIR, k, j, i) = Jx;
      J(JDIR, k, j, i) = Jy;
      J(KDIR, k, j, i) = Jz;
    }
  );

  // Regularize the current along the axis in cases where an axis is present.
  // This feature is not yet ready, so it is commented out
  if(this->haveAxis) {
    boundary->axis->RegularizeCurrent();
  }

  idfx::popRegion();
}
#endif //FLUID_CALCCURRENT_HPP_
