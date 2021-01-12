// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "../idefix.hpp"
#include "hydro.hpp"


// Build a left and right extrapolation of the primitive variables along direction dir

// These functions extrapolate the cell prim vars to the faces. Definitions are as followed
//
// |       cell i-1               interface i          cell i
// |-----------------------------------|------------------------------------||
//          Vc(i-1)           PrimL(i)  PrimR(i)       Vc(i)

void Hydro::ExtrapolatePrimVar(int dir) {
  int ioffset,joffset,koffset;
  int iextend, jextend,kextend;
  int BXn;

  idfx::pushRegion("Hydro::ExtrapolatePrimVar");
  // Offset is in the direction of integration
  ioffset=joffset=koffset=0;

  // extension if perp to the direction of integration, as required by CT.
  iextend=jextend=kextend=0;

  // Determine the offset along which we do the extrapolation, as well as the perp extension
  if(dir==IDIR) {
    ioffset=1;
    BXn = BX1;
    D_EXPAND(               ,
              jextend = 1;  ,
              kextend = 1;  )
  }
  if(dir==JDIR) {
    joffset=1;
    BXn = BX2;
    D_EXPAND( iextend = 1;  ,
                            ,
              kextend = 1;  )
  }
  if(dir==KDIR) {
    koffset=1;
    BXn = BX3;
    D_EXPAND( iextend = 1;  ,
              jextend = 1;  ,
                            )
  }

  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Vs = this->Vs;
  IdefixArray4D<real> PrimL = this->PrimL;
  IdefixArray4D<real> PrimR = this->PrimR;


#if ORDER == 1

  idefix_for("ExtrapolatePrimVar",
             0,NVAR,data->beg[KDIR]-kextend,data->end[KDIR]+koffset+kextend,
             data->beg[JDIR]-jextend,data->end[JDIR]+joffset+jextend,
             data->beg[IDIR]-iextend,data->end[IDIR]+ioffset+iextend,
    KOKKOS_LAMBDA (int n, int k, int j, int i) {
      // If normal component, the use Staggered field
      if(n==BXn) {
        PrimL(n,k,j,i) = Vs(dir,k,j,i);
        PrimR(n,k,j,i) = Vs(dir,k,j,i);
      } else {
        PrimL(n,k,j,i) = Vc(n,k-koffset,j-joffset,i-ioffset);
        PrimR(n,k,j,i) = Vc(n,k,j,i);
      }
  });

#elif ORDER == 2
  idefix_for("ExtrapolatePrimVar",
             0,NVAR,data->beg[KDIR]-koffset-kextend,data->end[KDIR]+koffset+kextend,
             data->beg[JDIR]-joffset-jextend,data->end[JDIR]+joffset+jextend,
             data->beg[IDIR]-ioffset-iextend,data->end[IDIR]+ioffset+iextend,
    KOKKOS_LAMBDA (int n, int k, int j, int i) {
      if(n==BXn) {
        PrimL(n,k+koffset,j+joffset,i+ioffset) = Vs(dir,k+koffset,j+joffset,i+ioffset);
        PrimR(n,k,j,i) = Vs(dir,k,j,i);
      } else {
        real dvm = Vc(n,k,j,i)-Vc(n,k-koffset,j-joffset,i-ioffset);
        real dvp = Vc(n,k+koffset,j+joffset,i+ioffset) - Vc(n,k,j,i);

        // Van Leer limiter
        real dv = (dvp*dvm > ZERO_F ? TWO_F*dvp*dvm/(dvp + dvm) : ZERO_F);

        PrimL(n,k+koffset,j+joffset,i+ioffset) = Vc(n,k,j,i) + HALF_F*dv;
        PrimR(n,k,j,i) = Vc(n,k,j,i) - HALF_F*dv;
      }
  });
#else
  #error ORDER should be 1 or 2
#endif


  idfx::popRegion();
}
