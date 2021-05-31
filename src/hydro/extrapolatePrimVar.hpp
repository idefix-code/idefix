// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************
#ifndef HYDRO_EXTRAPOLATEPRIMVAR_HPP_
#define HYDRO_EXTRAPOLATEPRIMVAR_HPP_

#include "hydro.hpp"
#include "dataBlock.hpp"

// Build a left and right extrapolation of the primitive variables along direction dir

// These functions extrapolate the cell prim vars to the faces. Definitions are as followed
//
// |       cell i-1               interface i          cell i
// |-----------------------------------|------------------------------------||
//          Vc(i-1)           PrimL(i)  PrimR(i)       Vc(i)

template<const int DIR>
KOKKOS_INLINE_FUNCTION void Hydro::K_ExtrapolatePrimVar
      (const int i, const int j, const int k, const IdefixArray4D<real> &Vc,
      const IdefixArray4D<real> &Vs, real vL[], real vR[]) {
  // 1-- Store the primitive variables on the left, right, and averaged states
  const int ioffset = (DIR==IDIR ? 1 : 0);
  const int joffset = (DIR==JDIR ? 1 : 0);
  const int koffset = (DIR==KDIR ? 1 : 0);
  const int BXn = BX1 + DIR;

  #pragma unroll
  for(int nv = 0 ; nv < NVAR; nv++) {
    #if ORDER == 1
      if(nv==BXn) {
        vR[nv] = vL[nv] = Vs(DIR,k,j,i);
      } else {
        vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset);
        vR[nv] = Vc(nv,k,j,i);
      }
    #elif ORDER == 2
      if(nv==BXn) {
        vR[nv] = vL[nv] = Vs(DIR,k,j,i);
      } else {
        real dvm = Vc(nv,k-koffset,j-joffset,i-ioffset)
                  -Vc(nv,k-2*koffset,j-2*joffset,i-2*ioffset);
        real dvp = Vc(nv,k,j,i)-Vc(nv,k-koffset,j-joffset,i-ioffset);
        // Van Leer limiter
        real dv = (dvp*dvm > ZERO_F ? TWO_F*dvp*dvm/(dvp + dvm) : ZERO_F);

        vL[nv] = Vc(nv,k-koffset,j-joffset,i-ioffset) + HALF_F*dv;

        dvm = dvp;
        dvp = Vc(nv,k+koffset,j+joffset,i+ioffset) - Vc(nv,k,j,i);

        // Van Leer limiter
        dv = (dvp*dvm > ZERO_F ? TWO_F*dvp*dvm/(dvp + dvm) : ZERO_F);

        vR[nv] = Vc(nv,k,j,i) - HALF_F*dv;
      }
    #endif
  }
}

#endif // HYDRO_EXTRAPOLATEPRIMVAR_HPP_
