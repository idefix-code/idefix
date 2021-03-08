// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "hydro.hpp"
#include "dataBlock.hpp"


KOKKOS_INLINE_FUNCTION real MC_LIM2 (const real dp, const real dm) {
  real dc, scrh;

  if (dp*dm < ZERO_F) return(ZERO_F);

  dc   = HALF_F*(dp + dm);
  scrh = TWO_F*(std::fabs(dp) < std::fabs(dm) ? dp:dm);
  return (std::fabs(dc) < std::fabs(scrh) ? dc:scrh);
}


void ElectroMotiveForce::calcRiemannEmf() {
  idfx::pushRegion("ElectroMotiveForce::calcRiemannEmf");

#if EMF_AVERAGE == UCT_HLL

  IdefixArray3D<real> dvx_dx = this->dvx_dx;
  IdefixArray3D<real> dvy_dx = this->dvy_dx;
  IdefixArray3D<real> dvx_dy = this->dvx_dy;
  IdefixArray3D<real> dvy_dy = this->dvy_dy;
  #if DIMENSIONS == 3
  IdefixArray3D<real> dvx_dz = this->dvx_dz;
  IdefixArray3D<real> dvy_dz = this->dvy_dz;
  IdefixArray3D<real> dvz_dz = this->dvz_dz;
  IdefixArray3D<real> dvz_dx = this->dvz_dx;
  IdefixArray3D<real> dvz_dy = this->dvz_dy;
  #endif

  IdefixArray3D<real> dbx_dy = this->dbx_dy;
  IdefixArray3D<real> dby_dx = this->dby_dx;
  #if DIMENSIONS == 3
  IdefixArray3D<real> dbz_dx = this->dbz_dx;
  IdefixArray3D<real> dbz_dy = this->dbz_dy;
  IdefixArray3D<real> dbx_dz = this->dbx_dz;
  IdefixArray3D<real> dby_dz = this->dby_dz;
  #endif

  IdefixArray3D<real> SxR = this->SxR;
  IdefixArray3D<real> SxL = this->SxL;
  IdefixArray3D<real> SyR = this->SyR;
  IdefixArray3D<real> SyL = this->SyL;
  IdefixArray3D<real> SzR = this->SzR;
  IdefixArray3D<real> SzL = this->SzL;

  // Corned EMFs
  IdefixArray3D<real> ex = this->ex;
  IdefixArray3D<real> ey = this->ey;
  IdefixArray3D<real> ez = this->ez;

  // Field
  IdefixArray4D<real> Vs = hydro->Vs;
  IdefixArray4D<real> Vc = hydro->Vc;

  idefix_for("CalcSlopesEMF",
             KOFFSET,data->np_tot[KDIR]-KOFFSET,
             JOFFSET,data->np_tot[JDIR]-JOFFSET,
             IOFFSET,data->np_tot[IDIR]-IOFFSET,
    KOKKOS_LAMBDA (int k, int j, int i) {
      dbx_dy(k,j,i) = MC_LIM2(Vs(BX1s,k,j+1,i) - Vs(BX1s,k,j,i),
                              Vs(BX1s,k,j,i)   - Vs(BX1s,k,j-1,i));
      dby_dx(k,j,i) = MC_LIM2(Vs(BX2s,k,j,i+1) - Vs(BX2s,k,j,i),
                              Vs(BX2s,k,j,i)   - Vs(BX2s,k,j,i-1));
    #if DIMENSIONS == 3
      dbx_dz(k,j,i) = MC_LIM2(Vs(BX1s,k+1,j,i) - Vs(BX1s,k,j,i),
                              Vs(BX1s,k,j,i)   - Vs(BX1s,k-1,j,i));
      dby_dz(k,j,i) = MC_LIM2(Vs(BX2s,k+1,j,i) - Vs(BX2s,k,j,i),
                              Vs(BX2s,k,j,i)   - Vs(BX2s,k-1,j,i));
      dbz_dx(k,j,i) = MC_LIM2(Vs(BX3s,k,j,i+1) - Vs(BX3s,k,j,i),
                              Vs(BX3s,k,j,i)   - Vs(BX3s,k,j,i-1));
      dbz_dy(k,j,i) = MC_LIM2(Vs(BX3s,k,j+1,i) - Vs(BX3s,k,j,i),
                              Vs(BX3s,k,j,i)   - Vs(BX3s,k,j-1,i));
    #endif
    }
  );

  idefix_for("CalcHLLSolver",
             data->beg[KDIR],data->end[KDIR]+KOFFSET,
             data->beg[JDIR],data->end[JDIR]+JOFFSET,
             data->beg[IDIR],data->end[IDIR]+IOFFSET,
    KOKKOS_LAMBDA (int k, int j, int i) {
      // EMF: Z component at (i-1/2, j-1/2, k)
      int im = i-1, jm = j-1, km = k-1;

      real a_xp = std::fmax(SxR(k,jm,i), SxR(k,j,i));
      real a_xm = std::fmax(SxL(k,jm,i), SxL(k,j,i));
      real a_yp = std::fmax(SyR(k,j,im), SyR(k,j,i));
      real a_ym = std::fmax(SyL(k,j,im), SyL(k,j,i));

      real bS = (Vs(BX1s,k,jm,i) + HALF_F*(dbx_dy(k,jm,i)));
      real bW = (Vs(BX2s,k,j,im) + HALF_F*(dby_dx(k,j,im)));
      real bN = (Vs(BX1s,k,j,i) - HALF_F*(dbx_dy(k,j,i)));
      real bE = (Vs(BX2s,k,j,i) - HALF_F*(dby_dx(k,j,i)));

      real eSW, eSE, eNE, eNW;

      eSW =   (Vc(VX2,k,jm,im) + HALF_F*(dvy_dx(k,jm,im) + dvy_dy(k,jm,im)))*bS
            - (Vc(VX1,k,jm,im) + HALF_F*(dvx_dx(k,jm,im) + dvx_dy(k,jm,im)))*bW;
      eSE =   (Vc(VX2,k,jm,i) - HALF_F*(dvy_dx(k,jm,i) - dvy_dy(k,jm,i)))*bS
            - (Vc(VX1,k,jm,i) - HALF_F*(dvx_dx(k,jm,i) - dvx_dy(k,jm,i)))*bE;
      eNE =   (Vc(VX2,k,j,i) - HALF_F*(dvy_dx(k,j,i) + dvy_dy(k,j,i)))*bN
            - (Vc(VX1,k,j,i) - HALF_F*(dvx_dx(k,j,i) + dvx_dy(k,j,i)))*bE;
      eNW =   (Vc(VX2,k,j,im) + HALF_F*(dvy_dx(k,j,im) - dvy_dy(k,j,im)))*bN
            - (Vc(VX1,k,j,im) + HALF_F*(dvx_dx(k,j,im) - dvx_dy(k,j,im)))*bW;

      ez(k,j,i) = a_xp*a_yp*eSW + a_xm*a_yp*eSE + a_xm*a_ym*eNE + a_xp*a_ym*eNW;
      ez(k,j,i) /= (a_xp + a_xm)*(a_yp + a_ym);
      ez(k,j,i) -= a_yp*a_ym*(bN - bS)/(a_yp + a_ym);
      ez(k,j,i) += a_xp*a_xm*(bE - bW)/(a_xp + a_xm);

  #if DIMENSIONS == 3
      // EMF: X component at (i, j-1/2, k-1/2)

      a_xp = std::fmax(SyR(km,j,i), SyR(k,j,i));
      a_xm = std::fmax(SyL(km,j,i), SyL(k,j,i));
      a_yp = std::fmax(SzR(k,jm,i), SzR(k,j,i));
      a_ym = std::fmax(SzL(k,jm,i), SzL(k,j,i));

      bS = (Vs(BX2s,km,j,i) + HALF_F*(dby_dz(km,j,i)));
      bW = (Vs(BX3s,k,jm,i) + HALF_F*(dbz_dy(k,jm,i)));
      bN = (Vs(BX2s,k,j,i) - HALF_F*(dby_dz(k,j,i)));
      bE = (Vs(BX3s,k,j,i) - HALF_F*(dbz_dy(k,j,i)));

      eSW =   (Vc(VX3,km,jm,i) + HALF_F*(dvz_dy(km,jm,i) + dvz_dz(km,jm,i)))*bS
            - (Vc(VX2,km,jm,i) + HALF_F*(dvy_dy(km,jm,i) + dvy_dz(km,jm,i)))*bW;
      eSE =   (Vc(VX3,km,j,i) - HALF_F*(dvz_dy(km,j,i) - dvz_dz(km,j,i)))*bS
            - (Vc(VX2,km,j,i) - HALF_F*(dvy_dy(km,j,i) - dvy_dz(km,j,i)))*bE;
      eNE =   (Vc(VX3,k,j,i) - HALF_F*(dvz_dy(k,j,i) + dvz_dz(k,j,i)))*bN
            - (Vc(VX2,k,j,i) - HALF_F*(dvy_dy(k,j,i) + dvy_dz(k,j,i)))*bE;
      eNW =   (Vc(VX3,k,jm,i) + HALF_F*(dvz_dy(k,jm,i) - dvz_dz(k,jm,i)))*bN
            - (Vc(VX2,k,jm,i) + HALF_F*(dvy_dy(k,jm,i) - dvy_dz(k,jm,i)))*bW;

      ex(k,j,i) = a_xp*a_yp*eSW + a_xm*a_yp*eSE + a_xm*a_ym*eNE + a_xp*a_ym*eNW;
      ex(k,j,i) /= (a_xp + a_xm)*(a_yp + a_ym);
      ex(k,j,i) -= a_yp*a_ym*(bN - bS)/(a_yp + a_ym);
      ex(k,j,i) += a_xp*a_xm*(bE - bW)/(a_xp + a_xm);

      // EMF: Y component at (i-1/2, j, k-1/2)

      a_xp = std::fmax(SzR(k,j,im), SzR(k,j,i));
      a_xm = std::fmax(SzL(k,j,im), SzL(k,j,i));
      a_yp = std::fmax(SxR(km,j,i), SxR(k,j,i));
      a_ym = std::fmax(SxL(km,j,i), SxL(k,j,i));

      bS = (Vs(BX3s,k,j,im) + HALF_F*(dbz_dx(k,j,im)));
      bW = (Vs(BX1s,km,j,i) + HALF_F*(dbx_dz(km,j,i)));
      bN = (Vs(BX3s,k,j,i) - HALF_F*(dbz_dx(k,j,i)));
      bE = (Vs(BX1s,k,j,i) - HALF_F*(dbx_dz(k,j,i)));

      eSW =   (Vc(VX1,km,j,im) + HALF_F*(dvx_dz(km,j,im) + dvx_dx(km,j,im)))*bS
            - (Vc(VX3,km,j,im) + HALF_F*(dvz_dz(km,j,im) + dvz_dx(km,j,im)))*bW;
      eSE =   (Vc(VX1,k,j,im) - HALF_F*(dvx_dz(k,j,im) - dvx_dx(k,j,im)))*bS
            - (Vc(VX3,k,j,im) - HALF_F*(dvz_dz(k,j,im) - dvz_dx(k,j,im)))*bE;
      eNE =   (Vc(VX1,k,j,i) - HALF_F*(dvx_dz(k,j,i) + dvx_dx(k,j,i)))*bN
            - (Vc(VX3,k,j,i) - HALF_F*(dvz_dz(k,j,i) + dvz_dx(k,j,i)))*bE;
      eNW =   (Vc(VX1,km,j,i) + HALF_F*(dvx_dz(km,j,i) - dvx_dx(km,j,i)))*bN
            - (Vc(VX3,km,j,i) + HALF_F*(dvz_dz(km,j,i) - dvz_dx(km,j,i)))*bW;

      ey(k,j,i) = a_xp*a_yp*eSW + a_xm*a_yp*eSE + a_xm*a_ym*eNE + a_xp*a_ym*eNW;
      ey(k,j,i) /= (a_xp + a_xm)*(a_yp + a_ym);
      ey(k,j,i) -= a_yp*a_ym*(bN - bS)/(a_yp + a_ym);
      ey(k,j,i) += a_xp*a_xm*(bE - bW)/(a_xp + a_xm);
  #endif  // DIMENSIONS == 3
    }
  );
#else

#endif  // EMF_AVERAGE

  idfx::popRegion();
}
