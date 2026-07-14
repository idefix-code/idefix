// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_CONSTRAINEDTRANSPORT_CALCRIEMANNEMF_HPP_
#define FLUID_CONSTRAINEDTRANSPORT_CALCRIEMANNEMF_HPP_

#include "fluid.hpp"
#include "dataBlock.hpp"
#include "slopeLimiter.hpp" // Contains the main slope limiters called in this function
#include "physics.hpp"

// A shortcut for the slopelimiter template parameters
// Note that we do not use the full class, so the template parameters are meaningless
// here
using SL = SlopeLimiter<>;

template<typename Phys>
void ConstrainedTransport<Phys>::CalcRiemannAverage() {
  idfx::pushRegion("ConstrainedTransport::calcRiemannAverage");

  // Corned EMFs
  IdefixArray3D<real> ez = this->ez;
#if DIMENSIONS == 3
  IdefixArray3D<real> ex = this->ex;
  IdefixArray3D<real> ey = this->ey;

  // Face-centered EMFs
  IdefixArray3D<real> exj = this->exj;
  IdefixArray3D<real> exk = this->exk;
  IdefixArray3D<real> eyi = this->eyi;
  IdefixArray3D<real> eyk = this->eyk;
#endif

  // Face-centered EMFs
  IdefixArray3D<real> ezi = this->ezi;
  IdefixArray3D<real> ezj = this->ezj;

  IdefixArray3D<real> axL = this->axL;
  IdefixArray3D<real> axR = this->axR;
  IdefixArray3D<real> ayL = this->ayL;
  IdefixArray3D<real> ayR = this->ayR;
#if DIMENSIONS == 3
  IdefixArray3D<real> azL = this->azL;
  IdefixArray3D<real> azR = this->azR;
#endif

  IdefixArray3D<real> dxL = this->dxL;
  IdefixArray3D<real> dxR = this->dxR;
  IdefixArray3D<real> dyL = this->dyL;
  IdefixArray3D<real> dyR = this->dyR;
#if DIMENSIONS == 3
  IdefixArray3D<real> dzL = this->dzL;
  IdefixArray3D<real> dzR = this->dzR;
#endif

  IdefixArray4D<real> Vs = hydro->Vs;

  idefix_for("CalcCenterEMF",
             data->beg[KDIR],data->end[KDIR]+KOFFSET,
             data->beg[JDIR],data->end[JDIR]+JOFFSET,
             data->beg[IDIR],data->end[IDIR]+IOFFSET,
    KOKKOS_LAMBDA (int k, int j, int i) {
      // IDIR
#if DIMENSIONS >= 2
      [[maybe_unused]] real phi, vL, vR, dv, bL, bR, db;
      [[maybe_unused]] real aL, aR, dL, dR;

      [[maybe_unused]] const int im = i-1, jm = j-1, km = k-1;

      // EMF: Z component at (i-1/2, j-1/2, k)
      aL = HALF_F*(axL(k,jm,i) + axL(k,jm+1,i));
      aR = HALF_F*(axR(k,jm,i) + axR(k,jm+1,i));
      dL = HALF_F*(dxL(k,jm,i) + dxL(k,jm+1,i));
      dR = HALF_F*(dxR(k,jm,i) + dxR(k,jm+1,i));

      #if ORDER == 4
        SL::getPPMStates( Vs(BX2s,k,j,im-2),
                          Vs(BX2s,k,j,im-1),
                          Vs(BX2s,k,j,im),
                          Vs(BX2s,k,j,im+1),
                          Vs(BX2s,k,j,im+2),
                          db, bL); // We don't use the left value at m1

        SL::getPPMStates( Vs(BX2s,k,j,i-2),
                          Vs(BX2s,k,j,i-1),
                          Vs(BX2s,k,j,i),
                          Vs(BX2s,k,j,i+1),
                          Vs(BX2s,k,j,i+2),
                          bR, db); // We don't use the right value at i
      #else
        db = SL::McLim(Vs(BX2s,k,j,im+1) - Vs(BX2s,k,j,im),
                       Vs(BX2s,k,j,im)   - Vs(BX2s,k,j,im-1));
        bL = Vs(BX2s,k,j,im) + HALF_F*db;

        db = SL::McLim(Vs(BX2s,k,j,i+1) - Vs(BX2s,k,j,i),
                       Vs(BX2s,k,j,i)   - Vs(BX2s,k,j,i-1));
        bR = Vs(BX2s,k,j,i) - HALF_F*db;

      #endif

      #if ORDER == 4
        SL::getPPMStates( ezj(k,j,im-2),
                          ezj(k,j,im-1),
                          ezj(k,j,im),
                          ezj(k,j,im+1),
                          ezj(k,j,im+2),
                          dv, vL); // We don't use the left value at m1

        SL::getPPMStates( ezj(k,j,i-2),
                          ezj(k,j,i-1),
                          ezj(k,j,i),
                          ezj(k,j,i+1),
                          ezj(k,j,i+2),
                          vR, dv); // We don't use the right value at i
      #else
        dv = SL::McLim(ezj(k,j,im+1) - ezj(k,j,im),
                       ezj(k,j,im)   - ezj(k,j,im-1));
        vL = ezj(k,j,im) + HALF_F*dv;

        dv = SL::McLim(ezj(k,j,i+1) - ezj(k,j,i),
                       ezj(k,j,i)   - ezj(k,j,i-1));
        vR = ezj(k,j,i) - HALF_F*dv;
      #endif

      phi = dR*bR - dL*bL;
      ez(k,j,i) = (aL*vL*bL + aR*vR*bR) + phi;
#endif

#if DIMENSIONS == 3
      // EMF: Y component at (i-1/2, j, k-1/2)
      aL = HALF_F*(axL(km,j,i) + axL(km+1,j,i));
      aR = HALF_F*(axR(km,j,i) + axR(km+1,j,i));
      dL = HALF_F*(dxL(km,j,i) + dxL(km+1,j,i));
      dR = HALF_F*(dxR(km,j,i) + dxR(km+1,j,i));

      #if ORDER == 4
        SL::getPPMStates( Vs(BX3s,k,j,im-2),
                          Vs(BX3s,k,j,im-1),
                          Vs(BX3s,k,j,im),
                          Vs(BX3s,k,j,im+1),
                          Vs(BX3s,k,j,im+2),
                          db, bL); // We don't use the left value at m1

        SL::getPPMStates( Vs(BX3s,k,j,i-2),
                          Vs(BX3s,k,j,i-1),
                          Vs(BX3s,k,j,i),
                          Vs(BX3s,k,j,i+1),
                          Vs(BX3s,k,j,i+2),
                          bR, db); // We don't use the right value at i
      #else
        db = SL::McLim(Vs(BX3s,k,j,im+1) - Vs(BX3s,k,j,im),
                       Vs(BX3s,k,j,im)   - Vs(BX3s,k,j,im-1));
        bL = Vs(BX3s,k,j,im) + HALF_F*db;

        db = SL::McLim(Vs(BX3s,k,j,i+1) - Vs(BX3s,k,j,i),
                       Vs(BX3s,k,j,i)   - Vs(BX3s,k,j,i-1));
        bR = Vs(BX3s,k,j,i) - HALF_F*db;
      #endif

      #if ORDER == 4
        SL::getPPMStates( eyk(k,j,im-2),
                          eyk(k,j,im-1),
                          eyk(k,j,im),
                          eyk(k,j,im+1),
                          eyk(k,j,im+2),
                          dv, vL); // We don't use the left value at m1

        SL::getPPMStates( eyk(k,j,i-2),
                          eyk(k,j,i-1),
                          eyk(k,j,i),
                          eyk(k,j,i+1),
                          eyk(k,j,i+2),
                          vR, dv); // We don't use the right value at i
      #else
        dv = SL::McLim(eyk(k,j,im+1) - eyk(k,j,im),
                       eyk(k,j,im)   - eyk(k,j,im-1));
        vL = eyk(k,j,im) + HALF_F*dv;

        dv = SL::McLim(eyk(k,j,i+1) - eyk(k,j,i),
                       eyk(k,j,i)   - eyk(k,j,i-1));
        vR = eyk(k,j,i) - HALF_F*dv;
      #endif

      phi = dR*bR - dL*bL;
      ey(k,j,i) = (aL*vL*bL + aR*vR*bR) - phi;
#endif

      // JDIR
#if DIMENSIONS >= 2
  #if DIMENSIONS == 3
      // EMF: X component at (i, j-1/2, k-1/2)
      aL = HALF_F*(ayL(km,j,i) + ayL(km+1,j,i));
      aR = HALF_F*(ayR(km,j,i) + ayR(km+1,j,i));
      dL = HALF_F*(dyL(km,j,i) + dyL(km+1,j,i));
      dR = HALF_F*(dyR(km,j,i) + dyR(km+1,j,i));

      #if ORDER == 4
        SL::getPPMStates(Vs(BX3s,k,jm-2,i),
                     Vs(BX3s,k,jm-1,i),
                     Vs(BX3s,k,jm,i),
                     Vs(BX3s,k,jm+1,i),
                     Vs(BX3s,k,jm+2,i),
                     db, bL);

        SL::getPPMStates(Vs(BX3s,k,j-2,i),
                     Vs(BX3s,k,j-1,i),
                     Vs(BX3s,k,j,i),
                     Vs(BX3s,k,j+1,i),
                     Vs(BX3s,k,j+2,i),
                     bR, db);
      #else
        db = SL::McLim(Vs(BX3s,k,jm+1,i) - Vs(BX3s,k,jm,i),
                    Vs(BX3s,k,jm,i)   - Vs(BX3s,k,jm-1,i));
        bL = Vs(BX3s,k,jm,i) + HALF_F*db;

        db = SL::McLim(Vs(BX3s,k,j+1,i) - Vs(BX3s,k,j,i),
                    Vs(BX3s,k,j,i)   - Vs(BX3s,k,j-1,i));
        bR = Vs(BX3s,k,j,i) - HALF_F*db;
      #endif

      #if ORDER == 4
        SL::getPPMStates(exk(k,jm-2,i),
                     exk(k,jm-1,i),
                     exk(k,jm,i),
                     exk(k,jm+1,i),
                     exk(k,jm+2,i),
                     dv, vL);

        SL::getPPMStates(exk(k,j-2,i),
                     exk(k,j-1,i),
                     exk(k,j,i),
                     exk(k,j+1,i),
                     exk(k,j+2,i),
                     vR, dv);
      #else
        dv = SL::McLim(exk(k,jm+1,i) - exk(k,jm,i),
                    exk(k,jm,i)   - exk(k,jm-1,i));
        vL = exk(k,jm,i) + HALF_F*dv;

        dv = SL::McLim(exk(k,j+1,i) - exk(k,j,i),
                    exk(k,j,i)   - exk(k,j-1,i));
        vR = exk(k,j,i) - HALF_F*dv;
      #endif

      phi = dR*bR - dL*bL;
      ex(k,j,i) = (aL*vL*bL + aR*vR*bR) + phi;
  #endif

      // EMF: Z component at (i-1/2, j-1/2, k)
      aL = HALF_F*(ayL(k,j,im) + ayL(k,j,im+1));
      aR = HALF_F*(ayR(k,j,im) + ayR(k,j,im+1));
      dL = HALF_F*(dyL(k,j,im) + dyL(k,j,im+1));
      dR = HALF_F*(dyR(k,j,im) + dyR(k,j,im+1));

      #if ORDER == 4
        SL::getPPMStates(Vs(BX1s,k,jm-2,i),
                     Vs(BX1s,k,jm-1,i),
                     Vs(BX1s,k,jm,i),
                     Vs(BX1s,k,jm+1,i),
                     Vs(BX1s,k,jm+2,i),
                     db, bL);

        SL::getPPMStates(Vs(BX1s,k,j-2,i),
                     Vs(BX1s,k,j-1,i),
                     Vs(BX1s,k,j,i),
                     Vs(BX1s,k,j+1,i),
                     Vs(BX1s,k,j+2,i),
                     bR, db);
      #else
        db = SL::McLim(Vs(BX1s,k,jm+1,i) - Vs(BX1s,k,jm,i),
                    Vs(BX1s,k,jm,i)   - Vs(BX1s,k,jm-1,i));
        bL = Vs(BX1s,k,jm,i) + HALF_F*db;

        db = SL::McLim(Vs(BX1s,k,j+1,i) - Vs(BX1s,k,j,i),
                    Vs(BX1s,k,j,i)   - Vs(BX1s,k,j-1,i));
        bR = Vs(BX1s,k,j,i) - HALF_F*db;
      #endif

      #if ORDER == 4
        SL::getPPMStates(ezi(k,jm-2,i),
                     ezi(k,jm-1,i),
                     ezi(k,jm,i),
                     ezi(k,jm+1,i),
                     ezi(k,jm+2,i),
                     dv, vL);

        SL::getPPMStates(ezi(k,j-2,i),
                     ezi(k,j-1,i),
                     ezi(k,j,i),
                     ezi(k,j+1,i),
                     ezi(k,j+2,i),
                     vR, dv);
      #else
        dv = SL::McLim(ezi(k,jm+1,i) - ezi(k,jm,i),
                    ezi(k,jm,i)   - ezi(k,jm-1,i));
        vL = ezi(k,jm,i) + HALF_F*dv;

        dv = SL::McLim(ezi(k,j+1,i) - ezi(k,j,i),
                    ezi(k,j,i)   - ezi(k,j-1,i));
        vR = ezi(k,j,i) - HALF_F*dv;
      #endif

      phi = dR*bR - dL*bL;
      ez(k,j,i) += (aL*vL*bL + aR*vR*bR) - phi;
#endif

      // KDIR
#if DIMENSIONS == 3
      // EMF: Y component at (i-1/2, j, k-1/2)
      aL = HALF_F*(azL(k,j,im) + azL(k,j,im+1));
      aR = HALF_F*(azR(k,j,im) + azR(k,j,im+1));
      dL = HALF_F*(dzL(k,j,im) + dzL(k,j,im+1));
      dR = HALF_F*(dzR(k,j,im) + dzR(k,j,im+1));

      #if ORDER == 4
        SL::getPPMStates(Vs(BX1s,km-2,j,i),
                     Vs(BX1s,km-1,j,i),
                     Vs(BX1s,km,j,i),
                     Vs(BX1s,km+1,j,i),
                     Vs(BX1s,km+2,j,i),
                     db, bL);

        SL::getPPMStates(Vs(BX1s,k-2,j,i),
                     Vs(BX1s,k-1,j,i),
                     Vs(BX1s,k,j,i),
                     Vs(BX1s,k+1,j,i),
                     Vs(BX1s,k+2,j,i),
                     bR, db);
      #else
        db = SL::McLim(Vs(BX1s,km+1,j,i) - Vs(BX1s,km,j,i),
                      Vs(BX1s,km,j,i)   - Vs(BX1s,km-1,j,i));
        bL = Vs(BX1s,km,j,i) + HALF_F*db;

        db = SL::McLim(Vs(BX1s,k+1,j,i) - Vs(BX1s,k,j,i),
                      Vs(BX1s,k,j,i)   - Vs(BX1s,k-1,j,i));
        bR = Vs(BX1s,k,j,i) - HALF_F*db;
      #endif

      #if ORDER == 4
        SL::getPPMStates(eyi(km-2,j,i),
                     eyi(km-1,j,i),
                     eyi(km,j,i),
                     eyi(km+1,j,i),
                     eyi(km+2,j,i),
                     dv, vL);

        SL::getPPMStates(eyi(k-2,j,i),
                     eyi(k-1,j,i),
                     eyi(k,j,i),
                     eyi(k+1,j,i),
                     eyi(k+2,j,i),
                     vR, dv);
      #else
        dv = SL::McLim(eyi(km+1,j,i) - eyi(km,j,i),
                    eyi(km,j,i)   - eyi(km-1,j,i));
        vL = eyi(km,j,i) + HALF_F*dv;

        dv = SL::McLim(eyi(k+1,j,i) - eyi(k,j,i),
                    eyi(k,j,i)   - eyi(k-1,j,i));
        vR = eyi(k,j,i) - HALF_F*dv;
      #endif

      phi = dR*bR - dL*bL;
      ey(k,j,i) += (aL*vL*bL + aR*vR*bR) + phi;

      // EMF: X component at (i, j-1/2, k-1/2)
      aL = HALF_F*(azL(k,jm,i) + azL(k,jm+1,i));
      aR = HALF_F*(azR(k,jm,i) + azR(k,jm+1,i));
      dL = HALF_F*(dzL(k,jm,i) + dzL(k,jm+1,i));
      dR = HALF_F*(dzR(k,jm,i) + dzR(k,jm+1,i));

      #if ORDER == 4
        SL::getPPMStates(Vs(BX2s,km-2,j,i),
                     Vs(BX2s,km-1,j,i),
                     Vs(BX2s,km,j,i),
                     Vs(BX2s,km+1,j,i),
                     Vs(BX2s,km+2,j,i),
                     db, bL);

        SL::getPPMStates(Vs(BX2s,k-2,j,i),
                     Vs(BX2s,k-1,j,i),
                     Vs(BX2s,k,j,i),
                     Vs(BX2s,k+1,j,i),
                     Vs(BX2s,k+2,j,i),
                     bR, db);
      #else
        db = SL::McLim(Vs(BX2s,km+1,j,i) - Vs(BX2s,km,j,i),
                    Vs(BX2s,km,j,i)   - Vs(BX2s,km-1,j,i));
        bL = Vs(BX2s,km,j,i) + HALF_F*db;

        db = SL::McLim(Vs(BX2s,k+1,j,i) - Vs(BX2s,k,j,i),
                    Vs(BX2s,k,j,i)   - Vs(BX2s,k-1,j,i));
        bR = Vs(BX2s,k,j,i) - HALF_F*db;
      #endif

      #if ORDER == 4
        SL::getPPMStates(exj(km-2,j,i),
                     exj(km-1,j,i),
                     exj(km,j,i),
                     exj(km+1,j,i),
                     exj(km+2,j,i),
                     dv, vL);

        SL::getPPMStates(exj(k-2,j,i),
                     exj(k-1,j,i),
                     exj(k,j,i),
                     exj(k+1,j,i),
                     exj(k+2,j,i),
                     vR, dv);
      #else
        dv = SL::McLim(exj(km+1,j,i) - exj(km,j,i),
                    exj(km,j,i)   - exj(km-1,j,i));
        vL = exj(km,j,i) + HALF_F*dv;

        dv = SL::McLim(exj(k+1,j,i) - exj(k,j,i),
                    exj(k,j,i)   - exj(k-1,j,i));
        vR = exj(k,j,i) - HALF_F*dv;
      #endif

      phi = dR*bR - dL*bL;
      ex(k,j,i) += (aL*vL*bL + aR*vR*bR) - phi;
#endif
    }
  );

  idfx::popRegion();
}
#endif // FLUID_CONSTRAINEDTRANSPORT_CALCRIEMANNEMF_HPP_
