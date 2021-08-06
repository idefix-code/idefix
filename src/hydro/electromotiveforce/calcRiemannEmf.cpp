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
/*
KOKKOS_INLINE_FUNCTION real LIMITER (const real dp, real dm) {

}
*/

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

    dvx_dx(k,j,i) = MC_LIM2(Vc(VX1,k,j,i+1) - Vc(VX1,k,j,i),
                            Vc(VX1,k,j,i)   - Vc(VX1,k,j,i-1));
    #if DIMENSIONS >= 2
      dvx_dy(k,j,i) = MC_LIM2(Vc(VX1,k,j+1,i) - Vc(VX1,k,j,i),
                              Vc(VX1,k,j,i)   - Vc(VX1,k,j-1,i));
    #endif
    #if COMPONENTS >= 2
      dvy_dx(k,j,i) = MC_LIM2(Vc(VX2,k,j,i+1) - Vc(VX2,k,j,i),
                              Vc(VX2,k,j,i)   - Vc(VX2,k,j,i-1));
      #if DIMENSIONS >= 2
      dvy_dy(k,j,i) = MC_LIM2(Vc(VX2,k,j+1,i) - Vc(VX2,k,j,i),
                              Vc(VX2,k,j,i)   - Vc(VX2,k,j-1,i));
      #endif
    #endif // COMPONENTS
    #if COMPONENTS == 3
      dvz_dx(k,j,i) = MC_LIM2(Vc(VX3,k,j,i+1) - Vc(VX3,k,j,i),
                              Vc(VX3,k,j,i)   - Vc(VX3,k,j,i-1));
      #if DIMENSIONS >= 2
      dvz_dy(k,j,i) = MC_LIM2(Vc(VX3,k,j+1,i) - Vc(VX3,k,j,i),
                              Vc(VX3,k,j,i)   - Vc(VX3,k,j-1,i));
      #endif
    #endif // COMPONENTS
    #if DIMENSIONS == 3
      dvx_dz(k,j,i) = MC_LIM2(Vc(VX1,k+1,j,i) - Vc(VX1,k,j,i),
                            Vc(VX1,k,j,i)   - Vc(VX1,k-1,j,i));
      dvy_dz(k,j,i) = MC_LIM2(Vc(VX2,k+1,j,i) - Vc(VX2,k,j,i),
                            Vc(VX2,k,j,i)   - Vc(VX2,k-1,j,i));
      dvz_dz(k,j,i) = MC_LIM2(Vc(VX3,k+1,j,i) - Vc(VX3,k,j,i),
                            Vc(VX3,k,j,i)   - Vc(VX3,k-1,j,i));
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

void ElectroMotiveForce::calcRiemann2DEmf() {
  idfx::pushRegion("ElectroMotiveForce::calcRiemann2DEmf");

  // Corned EMFs
  IdefixArray3D<real> ex = this->ex;
  IdefixArray3D<real> ey = this->ey;
  IdefixArray3D<real> ez = this->ez;

  // Face-centered EMFs
  IdefixArray3D<real> exj = this->exj;
  IdefixArray3D<real> exk = this->exk;
  IdefixArray3D<real> eyi = this->eyi;
  IdefixArray3D<real> eyk = this->eyk;
  IdefixArray3D<real> ezi = this->ezi;
  IdefixArray3D<real> ezj = this->ezj;

  IdefixArray3D<real> axL = this->axL;
  IdefixArray3D<real> axR = this->axR;
  IdefixArray3D<real> ayL = this->ayL;
  IdefixArray3D<real> ayR = this->ayR;
  IdefixArray3D<real> azL = this->azL;
  IdefixArray3D<real> azR = this->azR;

  IdefixArray3D<real> dxL = this->dxL;
  IdefixArray3D<real> dxR = this->dxR;
  IdefixArray3D<real> dyL = this->dyL;
  IdefixArray3D<real> dyR = this->dyR;
  IdefixArray3D<real> dzL = this->dzL;
  IdefixArray3D<real> dzR = this->dzR;

  IdefixArray4D<real> Vs = hydro->Vs;

  idefix_for("CalcCenterEMF",
             2*KOFFSET,data->np_tot[KDIR]-2*KOFFSET,
             2*JOFFSET,data->np_tot[JDIR]-2*JOFFSET,
             2*IOFFSET,data->np_tot[IDIR]-2*IOFFSET,
    KOKKOS_LAMBDA (int k, int j, int i) {
      real phi, vL, vR, dv, bL, bR, db;
      real aL, aR, dL, dR;

      // IDIR
#if DIMENSIONS >= 2
      aL = HALF_F*(axL(k,j-1,i) + axL(k,j,i));
      aR = HALF_F*(axR(k,j-1,i) + axR(k,j,i));
      dL = HALF_F*(dxL(k,j-1,i) + dxL(k,j,i));
      dR = HALF_F*(dxR(k,j-1,i) + dxR(k,j,i));

      db = MC_LIM2(Vs(BX2s,k,j,i) - Vs(BX2s,k,j,i-1),
                   Vs(BX2s,k,j,i-1)   - Vs(BX2s,k,j,i-2));
      bL = Vs(BX2s,k,j,i-1) + HALF_F*db;

      db = MC_LIM2(Vs(BX2s,k,j,i+1) - Vs(BX2s,k,j,i),
                   Vs(BX2s,k,j,i)   - Vs(BX2s,k,j,i-1));
      bR = Vs(BX2s,k,j,i) - HALF_F*db;

      dv = MC_LIM2(ezj(k,j,i) - ezj(k,j,i-1),
                   ezj(k,j,i-1)   - ezj(k,j,i-2));
      vL = ezj(k,j,i-1) + HALF_F*dv;

      dv = MC_LIM2(ezj(k,j,i+1) - ezj(k,j,i),
                   ezj(k,j,i)   - ezj(k,j,i-1));
      vR = ezj(k,j,i) - HALF_F*dv;

      phi = dR*bR - dL*bL;
      ez(k,j,i) = (aL*vL*bL + aR*vR*bR) + phi;
#endif

#if DIMENSIONS == 3
      aL = HALF_F*(axL(k-1,j,i) + axL(k,j,i));
      aR = HALF_F*(axR(k-1,j,i) + axR(k,j,i));
      dL = HALF_F*(dxL(k-1,j,i) + dxL(k,j,i));
      dR = HALF_F*(dxR(k-1,j,i) + dxR(k,j,i));

      db = MC_LIM2(Vs(BX3s,k,j,i+1) - Vs(BX3s,k,j,i),
                   Vs(BX3s,k,j,i)   - Vs(BX3s,k,j,i-1));
      bL = Vs(BX3s,k,j,i) + HALF_F*db;
      bR = Vs(BX3s,k,j,i) - HALF_F*db;

      dv = MC_LIM2(eyk(k,j,i+1) - eyk(k,j,i),
                   eyk(k,j,i)   - eyk(k,j,i-1));
      vL = eyk(k,j,i) + HALF_F*dv;
      vR = eyk(k,j,i) - HALF_F*dv;

      phi = dR*bR - dL*bL;
      ey(k,j,i) = (aL*vL*bL + aR*vR*bR) - phi;
#endif

      // JDIR
#if DIMENSIONS >= 2
  #if DIMENSIONS == 3
      aL = HALF_F*(ayL(k-1,j,i) + ayL(k,j,i));
      aR = HALF_F*(ayR(k-1,j,i) + ayR(k,j,i));
      dL = HALF_F*(dyL(k-1,j,i) + dyL(k,j,i));
      dR = HALF_F*(dyR(k-1,j,i) + dyR(k,j,i));

      db = MC_LIM2(Vs(BX3s,k,j+1,i) - Vs(BX3s,k,j,i),
                   Vs(BX3s,k,j,i)   - Vs(BX3s,k,j-1,i));
      bL = Vs(BX3s,k,j,i) + HALF_F*db;
      bR = Vs(BX3s,k,j,i) - HALF_F*db;

      dv = MC_LIM2(exk(k,j+1,i) - exk(k,j,i),
                   exk(k,j,i)   - exk(k,j-1,i));
      vL = exk(k,j,i) + HALF_F*dv;
      vR = exk(k,j,i) - HALF_F*dv;

      phi = dR*bR - dL*bL;
      ex(k,j,i) = (aL*vL*bL + aR*vR*bR) + phi;
  #endif

      aL = HALF_F*(ayL(k,j,i-1) + ayL(k,j,i));
      aR = HALF_F*(ayR(k,j,i-1) + ayR(k,j,i));
      dL = HALF_F*(dyL(k,j,i-1) + dyL(k,j,i));
      dR = HALF_F*(dyR(k,j,i-1) + dyR(k,j,i));

      db = MC_LIM2(Vs(BX1s,k,j,i) - Vs(BX1s,k,j-1,i),
                   Vs(BX1s,k,j-1,i)  - Vs(BX1s,k,j-2,i));
      bL = Vs(BX1s,k,j-1,i) + HALF_F*db;

      db = MC_LIM2(Vs(BX1s,k,j+1,i) - Vs(BX1s,k,j,i),
                   Vs(BX1s,k,j,i)   - Vs(BX1s,k,j-1,i));
      bR = Vs(BX1s,k,j,i) - HALF_F*db;

      dv = MC_LIM2(ezi(k,j,i) - ezi(k,j-1,i),
                   ezi(k,j-1,i)   - ezi(k,j-2,i));
      vL = ezi(k,j-1,i) + HALF_F*dv;

      dv = MC_LIM2(ezi(k,j+1,i) - ezi(k,j,i),
                   ezi(k,j,i)   - ezi(k,j-1,i));
      vR = ezi(k,j,i) - HALF_F*dv;

      phi = dR*bR - dL*bL;
      ez(k,j,i) += (aL*vL*bL + aR*vR*bR) - phi;
#endif

      // KDIR
#if DIMENSIONS == 3
      aL = HALF_F*(azL(k,j,i-1) + azL(k,j,i));
      aR = HALF_F*(azR(k,j,i-1) + azR(k,j,i));
      dL = HALF_F*(dzL(k,j,i-1) + dzL(k,j,i));
      dR = HALF_F*(dzR(k,j,i-1) + dzR(k,j,i));

      db = MC_LIM2(Vs(BX1s,k+1,j,i) - Vs(BX1s,k,j,i),
                   Vs(BX1s,k,j,i)   - Vs(BX1s,k-1,j,i));
      bL = Vs(BX1s,k,j,i) + HALF_F*db;
      bR = Vs(BX1s,k,j,i) - HALF_F*db;

      dv = MC_LIM2(eyi(k+1,j,i) - eyi(k,j,i),
                   eyi(k,j,i)   - eyi(k-1,j,i));
      vL = eyi(k,j,i) + HALF_F*dv;
      vR = eyi(k,j,i) - HALF_F*dv;

      phi = dR*bR - dL*bL;
      ey(k,j,i) += (aL*vL*bL + aR*vR*bR) + phi;

  #if DIMENSIONS >= 2
      aL = HALF_F*(azL(k,j-1,i) + azL(k,j,i));
      aR = HALF_F*(azR(k,j-1,i) + azR(k,j,i));
      dL = HALF_F*(dzL(k,j-1,i) + dzL(k,j,i));
      dR = HALF_F*(dzR(k,j-1,i) + dzR(k,j,i));

      db = MC_LIM2(Vs(BX2s,k+1,j,i) - Vs(BX2s,k,j,i),
                   Vs(BX2s,k,j,i)   - Vs(BX2s,k-1,j,i));
      bL = Vs(BX2s,k,j,i) + HALF_F*db;
      bR = Vs(BX2s,k,j,i) - HALF_F*db;

      dv = MC_LIM2(exj(k+1,j,i) - exj(k,j,i),
                   exj(k,j,i)   - exj(k-1,j,i));
      vL = exj(k,j,i) + HALF_F*dv;
      vR = exj(k,j,i) - HALF_F*dv;

      phi = dR*bR - dL*bL;
      ex(k,j,i) += (aL*vL*bL + aR*vR*bR) - phi;
  #endif
#endif
    }
  );

  idfx::popRegion();
}
