// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_MHDSOLVERS_HLLDMHD_HPP_
#define HYDRO_MHDSOLVERS_HLLDMHD_HPP_

#include "../idefix.hpp"
#include "solversMHD.hpp"
#include "extrapolatePrimVar.hpp"
#include "fluxMHD.hpp"
#include "convertConsToPrimMHD.hpp"

// Compute Riemann fluxes from states using HLLD solver
template<const int DIR>
void Hydro::HlldMHD() {
  idfx::pushRegion("Hydro::HLLD_MHD");

  int ioffset,joffset,koffset;
  int iextend, jextend,kextend;
  ioffset=joffset=koffset=0;
  // extension in perp to the direction of integration, as required by CT.
  iextend=jextend=kextend=0;
  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Vs = this->Vs;
  IdefixArray4D<real> Flux = this->FluxRiemann;
  IdefixArray3D<real> cMax = this->cMax;
  IdefixArray3D<real> csIsoArr = this->isoSoundSpeedArray;

  // References to required emf components
  IdefixArray3D<real> Eb;
  IdefixArray3D<real> Et;

#if EMF_AVERAGE == UCT_CONTACT
  IdefixArray3D<int> SV;
#elif EMF_AVERAGE == UCT_HLL
  IdefixArray3D<real> SL;
  IdefixArray3D<real> SR;
#endif

  real gamma = this->gamma;
  real gamma_m1 = this->gamma-ONE_F;
  real csIso = this->isoSoundSpeed;
  HydroModuleStatus haveIsoCs = this->haveIsoSoundSpeed;


  // st and sb will be useful only when Hall is included
  D_EXPAND( real st;  ,
                      ,
            real sb;  )

  switch(DIR) {
    case(IDIR):
      ioffset = 1;
      D_EXPAND(
                st = -ONE_F;  ,
                jextend = 1;  ,
                kextend = 1;
                sb = +ONE_F;  )

      Et = this->emf.ezi;
      Eb = this->emf.eyi;
#if EMF_AVERAGE == UCT_CONTACT
      SV = this->emf.svx;
#elif EMF_AVERAGE == UCT_HLL
      SL = this->emf.SxL;
      SR = this->emf.SxR;
#endif

      break;
    case(JDIR):
      joffset=1;
      D_EXPAND(
                iextend = 1;
                st = +ONE_F;  ,
                              ,
                kextend = 1;
                sb = -ONE_F;  )

      Et = this->emf.ezj;
      Eb = this->emf.exj;
#if EMF_AVERAGE == UCT_CONTACT
      SV = this->emf.svy;
#elif EMF_AVERAGE == UCT_HLL
      SL = this->emf.SyL;
      SR = this->emf.SyR;
#endif

      break;
    case(KDIR):
      koffset=1;
      D_EXPAND(
                iextend = 1;
                st = -ONE_F;  ,
                jextend = 1;  ,
                sb = +ONE_F;  )

      Et = this->emf.eyk;
      Eb = this->emf.exk;
#if EMF_AVERAGE == UCT_CONTACT
      SV = this->emf.svz;
#elif EMF_AVERAGE == UCT_HLL
      SL = this->emf.SzL;
      SR = this->emf.SzR;
#endif

      break;
    default:
      IDEFIX_ERROR("Wrong direction");
  }

  idefix_for("CalcRiemannFlux",
             data->beg[KDIR]-kextend,data->end[KDIR]+koffset+kextend,
             data->beg[JDIR]-jextend,data->end[JDIR]+joffset+jextend,
             data->beg[IDIR]-iextend,data->end[IDIR]+ioffset+iextend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      // Init the directions (should be in the kernel for proper optimisation by the compilers)
      EXPAND( const int Xn = DIR+MX1;                    ,
              const int Xt = (DIR == IDIR ? MX2 : MX1);  ,
              const int Xb = (DIR == KDIR ? MX2 : MX3);  )

      EXPAND( const int BXn = DIR+BX1;                    ,
              const int BXt = (DIR == IDIR ? BX2 : BX1);  ,
              const int BXb = (DIR == KDIR ? BX2 : BX3);   )

      // Primitive variables
      real vL[NVAR];
      real vR[NVAR];

      K_ExtrapolatePrimVar<DIR>(i, j, k, Vc, Vs, vL, vR);

      // Conservative variables
      real uL[NVAR];
      real uR[NVAR];

      // Flux (left and right)
      real fluxL[NVAR];
      real fluxR[NVAR];

      // Signal speeds
      real cL, cR, cmax, c2Iso;

      // Init c2Isothermal (used only when isothermal approx is set)
      c2Iso = ZERO_F;

      // 2-- Get the wave speed
      real gpr, b1, b2, b3, Btmag2, Bmag2;
#if HAVE_ENERGY
      gpr = gamma*vL[PRS];
#else
      if(haveIsoCs == UserDefFunction) {
        c2Iso = HALF_F*(csIsoArr(k,j,i)+csIsoArr(k-koffset,j-joffset,i-ioffset));
        c2Iso = c2Iso*c2Iso;
      } else {
        c2Iso = csIso*csIso;
      }

      gpr = c2Iso*vL[RHO];
#endif

      // -- get total field
      b1 = b2 = b3 = ZERO_F;
      EXPAND ( b1 = vL[BXn];  ,
               b2 = vL[BXt];  ,
               b3 = vL[BXb];  )

      Btmag2 = b2*b2 + b3*b3;
      Bmag2  = b1*b1 + Btmag2;

      cL = gpr - Bmag2;
      cL = gpr + Bmag2 + std::sqrt(cL*cL + FOUR_F*gpr*Btmag2);
      cL = std::sqrt(HALF_F*cL/vL[RHO]);

#if HAVE_ENERGY
      gpr = gamma*vR[PRS];
#else
      gpr = c2Iso*vR[RHO];
#endif

      // -- get total field
      b1 = b2 = b3 = ZERO_F;
      EXPAND ( b1 = vR[BXn];  ,
               b2 = vR[BXt];  ,
               b3 = vR[BXb];  )

      Btmag2 = b2*b2 + b3*b3;
      Bmag2  = b1*b1 + Btmag2;

      cR = gpr - Bmag2;
      cR = gpr + Bmag2 + std::sqrt(cR*cR + FOUR_F*gpr*Btmag2);
      cR = std::sqrt(HALF_F*cR/vR[RHO]);

      // 4.1
      real cminL = vL[Xn] - cL;
      real cmaxL = vL[Xn] + cL;

      real cminR = vR[Xn] - cR;
      real cmaxR = vR[Xn] + cR;

      real sl = FMIN(cminL, cminR);
      real sr = FMAX(cmaxL, cmaxR);

      cmax  = std::fmax(FABS(sl), FABS(sr));

      // 2-- Compute the conservative variables
      K_PrimToCons(uL, vL, gamma_m1);
      K_PrimToCons(uR, vR, gamma_m1);

      // 3-- Compute the left and right fluxes
#pragma unroll
      for(int nv = 0 ; nv < NVAR; nv++) {
        fluxL[nv] = uL[nv];
        fluxR[nv] = uR[nv];
      }

      K_Flux(fluxL, vL, fluxL, c2Iso, ARG_EXPAND(Xn, Xt, Xb), ARG_EXPAND(BXn, BXt, BXb));
      K_Flux(fluxR, vR, fluxR, c2Iso, ARG_EXPAND(Xn, Xt, Xb), ARG_EXPAND(BXn, BXt, BXb));

      real ptR, ptL;

#if HAVE_ENERGY
      ptL  = vL[PRS] + 0.5* ( EXPAND(vL[BX1]*vL[BX1] , + vL[BX2]*vL[BX2], + vL[BX3]*vL[BX3]) );
      ptR  = vR[PRS] + 0.5* ( EXPAND(vR[BX1]*vR[BX1] , + vR[BX2]*vR[BX2], + vR[BX3]*vR[BX3]) );
#else
      ptL  = c2Iso*vL[RHO] + 0.5* (EXPAND(vL[BX1]*vL[BX1], + vL[BX2]*vL[BX2], + vL[BX3]*vL[BX3]));
      ptR  = c2Iso*vR[RHO] + 0.5* (EXPAND(vR[BX1]*vR[BX1], + vR[BX2]*vR[BX2], + vR[BX3]*vR[BX3]));
#endif

      // 5-- Compute the flux from the left and right states
      if (sl > 0) {
#pragma unroll
        for (int nv = 0 ; nv < NFLX; nv++) {
          Flux(nv,k,j,i) = fluxL[nv];
        }
      } else if (sr < 0) {
#pragma unroll
        for (int nv = 0 ; nv < NFLX; nv++) {
          Flux(nv,k,j,i) = fluxR[nv];
        }
      } else {
        real usL[NVAR];
        real usR[NVAR];

        real scrh, scrhL, scrhR, duL, duR, sBx, Bx, Bx1, SM, S1L, S1R;

#if HAVE_ENERGY
        real Uhll[NVAR];
        real vs, pts, sqrL, sqrR, vsL, vsR, wsL, wsR;
        int revert_to_hllc;

        // 3c. Compute U*(L), U^*(R)
        scrh = 1.0/(sr - sl);
        Bx1  = Bx = (sr*vR[BXn] - sl*vL[BXn])*scrh;
        sBx  = (Bx > 0.0 ? 1.0 : -1.0);

        duL  = sl - vL[Xn];
        duR  = sr - vR[Xn];

        scrh = 1.0/(duR*uR[RHO] - duL*uL[RHO]);
        SM   = (duR*uR[Xn] - duL*uL[Xn] - ptR + ptL)*scrh;

        pts  = duR*uR[RHO]*ptL - duL*uL[RHO]*ptR +
               vL[RHO]*vR[RHO]*duR*duL*(vR[Xn]- vL[Xn]);
        pts *= scrh;

        usL[RHO] = uL[RHO]*duL/(sl - SM);
        usR[RHO] = uR[RHO]*duR/(sr - SM);

        sqrL = std::sqrt(usL[RHO]);
        sqrR = std::sqrt(usR[RHO]);

        S1L = SM - fabs(Bx)/sqrL;
        S1R = SM + fabs(Bx)/sqrR;

        /* -----------------------------------------------------------------
        3d When S1L -> sl or S1R -> sr a degeneracy occurs.
        Although Miyoshi & Kusano say that no jump exists, we don't
        think this is actually true.
        Indeed, vy*, vz*, By*, Bz* cannot be solved independently.
        In this case we revert to the HLLC solver of Li (2005),  except
        for the term v.B in the region, which we compute in our own way.
        Note, that by comparing the expressions of Li (2005) and
        Miyoshi & Kusano (2005), the only change involves a
        re-definition of By* and Bz* in terms of By(HLL), Bz(HLL).
        ----------------------------------------------------------------- */

        revert_to_hllc = 0;

        if ( (S1L - sl) <  1.e-4*(SM - sl) ) revert_to_hllc = 1;
        if ( (S1R - sr) > -1.e-4*(sr - SM) ) revert_to_hllc = 1;

        if (revert_to_hllc) {
          scrh = 1.0/(sr - sl);
#pragma unroll
          for(int nv = 0 ; nv < NVAR; nv++) {
            Uhll[nv]  = sr*uR[nv] - sl*uL[nv] + fluxL[nv] - fluxR[nv];
            Uhll[nv] *= scrh;
          }

          // WHERE'S THE PRESSURE ?!?!?!?
          EXPAND( usL[BXn] = usR[BXn] = Uhll[BXn];  ,
                  usL[BXt] = usR[BXt] = Uhll[BXt];  ,
                  usL[BXb] = usR[BXb] = Uhll[BXb];  )

          S1L = S1R = SM; // region ** should never be computed since
                          // fluxes are given in terms of UL* and UR*
        } else {
          // 3e. Compute states in the * regions
          scrhL = (uL[RHO]*duL*duL - Bx*Bx)/(uL[RHO]*duL*(sl - SM) - Bx*Bx);
          scrhR = (uR[RHO]*duR*duR - Bx*Bx)/(uR[RHO]*duR*(sr - SM) - Bx*Bx);

          EXPAND( usL[BXn]  = Bx1;            ,
                  usL[BXt]  = uL[BXt]*scrhL;  ,
                  usL[BXb]  = uL[BXb]*scrhL;  )

          EXPAND( usR[BXn] = Bx1;            ,
                  usR[BXt] = uR[BXt]*scrhR;  ,
                  usR[BXb] = uR[BXb]*scrhR;  )
        }

        scrhL = Bx/(uL[RHO]*duL);
        scrhR = Bx/(uR[RHO]*duR);

        EXPAND(                                          ;  ,
                vsL = vL[Xt] - scrhL*(usL[BXt] - uL[BXt]);
                vsR = vR[Xt] - scrhR*(usR[BXt] - uR[BXt]);  ,

                wsL = vL[Xb] - scrhL*(usL[BXb] - uL[BXb]);
                wsR = vR[Xb] - scrhR*(usR[BXb] - uR[BXb]);  )

        EXPAND( usL[Xn] = usL[RHO]*SM;
                usR[Xn] = usR[RHO]*SM;   ,

                usL[Xt] = usL[RHO]*vsL;
                usR[Xt] = usR[RHO]*vsR;  ,

                usL[Xb] = usL[RHO]*wsL;
                usR[Xb] = usR[RHO]*wsR;  )

        /* -- Energy -- */

        scrhL  = EXPAND( vL[Xn]*Bx1, + vL[Xt]*uL[BXt], + vL[Xb]*uL[BXb]);
        scrhL -= EXPAND( SM*Bx1,     + vsL*usL[BXt],   + wsL*usL[BXb]);
        usL[ENG]  = duL*uL[ENG] - ptL*vL[Xn] + pts*SM + Bx*scrhL;
        usL[ENG] /= sl - SM;

        scrhR  = EXPAND(vR[Xn]*Bx1, + vR[Xt]*uR[BXt], + vR[Xb]*uR[BXb]);
        scrhR -= EXPAND(     SM*Bx1, +    vsR*usR[BXt], +    wsR*usR[BXb]);
        usR[ENG] = duR*uR[ENG] - ptR*vR[Xn] + pts*SM + Bx*scrhR;
        usR[ENG] /= sr - SM;

    // 3c. Compute flux when S1L > 0 or S1R < 0

        if (S1L >= 0.0) {       //  ----  Region L*
#pragma unroll
          for(int nv = 0 ; nv < NVAR; nv++) {
            Flux(nv,k,j,i) = fluxL[nv] + sl*(usL[nv] - uL[nv]);
          }
        } else if (S1R <= 0.0) {    //  ----  Region R*
#pragma unroll
          for(int nv = 0 ; nv < NVAR; nv++) {
            Flux(nv,k,j,i) = fluxR[nv] + sr*(usR[nv] - uR[nv]);
          }
        } else {   // -- This state exists only if B_x != 0
          // Compute U**
          real vss, wss;
          real ussl[NVAR];
          real ussr[NVAR];

          ussl[RHO] = usL[RHO];
          ussr[RHO] = usR[RHO];

          EXPAND(                      ,
                  vss  = sqrL*vsL + sqrR*vsR + (usR[BXt] - usL[BXt])*sBx;
                  vss /= sqrL + sqrR;  ,

                  wss  = sqrL*wsL + sqrR*wsR + (usR[BXb] - usL[BXb])*sBx;
                  wss /= sqrL + sqrR;  )

          EXPAND( ussl[Xn] = ussl[RHO]*SM;
                  ussr[Xn] = ussr[RHO]*SM;   ,

                  ussl[Xt] = ussl[RHO]*vss;
                  ussr[Xt] = ussr[RHO]*vss;  ,

                  ussl[Xb] = ussl[RHO]*wss;
                  ussr[Xb] = ussr[RHO]*wss;  )

          EXPAND( ussl[BXn] = ussr[BXn] = Bx1;  ,

                  ussl[BXt]  = sqrL*usR[BXt] + sqrR*usL[BXt] + sqrL*sqrR*(vsR - vsL)*sBx;
                  ussl[BXt] /= sqrL + sqrR;
                  ussr[BXt]  = ussl[BXt];       ,

                  ussl[BXb]  = sqrL*usR[BXb] + sqrR*usL[BXb] + sqrL*sqrR*(wsR - wsL)*sBx;
                  ussl[BXb] /= sqrL + sqrR;
                  ussr[BXb]  = ussl[BXb];      )

          // -- Energy jump

          scrhL  = EXPAND(SM*Bx1, +  vsL*usL[BXt], +  wsL*usL[BXb]);
          scrhL -= EXPAND(SM*Bx1, +  vss*ussl[BXt], +  wss*ussl[BXb]);

          scrhR  = EXPAND(SM*Bx1, +  vsR*usR[BXt], +  wsR*usR[BXb]);
          scrhR -= EXPAND(SM*Bx1, +  vss*ussr[BXt], +  wss*ussr[BXb]);

          ussl[ENG] = usL[ENG] - sqrL*scrhL*sBx;
          ussr[ENG] = usR[ENG] + sqrR*scrhR*sBx;


          if (SM >= 0.0) { //  ----  Region L**
#pragma unroll
            for(int nv = 0 ; nv < NVAR; nv++) {
              Flux(nv,k,j,i) = fluxL[nv] + S1L*(ussl[nv]  - usL[nv])
                              + sl*(usL[nv] - uL[nv]);
              }
          } else {         //  ----  Region R**
#pragma unroll
            for(int nv = 0 ; nv < NVAR; nv++) {
              Flux(nv,k,j,i) = fluxR[nv] + S1R*(ussr[nv]  - usR[nv])
                              + sr*(usR[nv] - uR[nv]);
            }
          }
        }  // end if (S1L < 0 S1R > 0)
#else
        real usc[NVAR];
        real rho, sqrho;
        int revert_to_hll;

        scrh = 1.0/(sr - sl);
        duL = sl - vL[Xn];
        duR = sr - vR[Xn];

        Bx1 = Bx = (sr*vR[BXn] - sl*vL[BXn])*scrh;

        rho                = (uR[RHO]*duR - uL[RHO]*duL)*scrh;
        Flux(RHO,k,j,i) = (sl*uR[RHO]*duR - sr*uL[RHO]*duL)*scrh;

        /* ---------------------------
            compute S*
        --------------------------- */

        sqrho = std::sqrt(rho);

        SM  = Flux(RHO,k,j,i)/rho;
        S1L = SM - fabs(Bx)/sqrho;
        S1R = SM + fabs(Bx)/sqrho;

        /* ---------------------------------------------
            Prevent degeneracies when S1L -> sl or
            S1R -> sr. Revert to HLL if necessary.
        --------------------------------------------- */

        revert_to_hll = 0;

        if ( (S1L - sl) <  1.e-4*(sr - sl) ) revert_to_hll = 1;
        if ( (S1R - sr) > -1.e-4*(sr - sl) ) revert_to_hll = 1;

        if (revert_to_hll) {
          scrh = 1.0/(sr - sl);
#pragma unroll
          for(int nv = 0 ; nv < NVAR; nv++) {
            Flux(nv,k,j,i) = sl*sr*(uR[nv] - uL[nv])
                            + sr*fluxL[nv] - sl*fluxR[nv];
            Flux(nv,k,j,i) *= scrh;
          }
        } else {
          Flux(Xn,k,j,i) = (sr*fluxL[Xn] - sl*fluxR[Xn]
                          + sr*sl*(uR[Xn] - uL[Xn]))*scrh;

          Flux(BXn,k,j,i) = sr*sl*(uR[BXn] - uL[BXn])*scrh;

      /* ---------------------------
                  Compute U*
          --------------------------- */

          scrhL = 1.0/((sl - S1L)*(sl - S1R));
          scrhR = 1.0/((sr - S1L)*(sr - S1R));

          EXPAND(                                                      ;  ,
                  usL[Xt] = rho*vL[Xt] - Bx*uL[BXt]*(SM - vL[Xn])*scrhL;
                  usR[Xt] = rho*vR[Xt] - Bx*uR[BXt]*(SM - vR[Xn])*scrhR;  ,

                  usL[Xb] = rho*vL[Xb] - Bx*uL[BXb]*(SM - vL[Xn])*scrhL;
                  usR[Xb] = rho*vR[Xb] - Bx*uR[BXb]*(SM - vR[Xn])*scrhR;  )

          EXPAND(                                                       ;  ,
                  usL[BXt] = uL[BXt]/rho*(uL[RHO]*duL*duL - Bx*Bx)*scrhL;
                  usR[BXt] = uR[BXt]/rho*(uR[RHO]*duR*duR - Bx*Bx)*scrhR;  ,

                  usL[BXb] = uL[BXb]/rho*(uL[RHO]*duL*duL - Bx*Bx)*scrhL;
                  usR[BXb] = uR[BXb]/rho*(uR[RHO]*duR*duR - Bx*Bx)*scrhR;  )

          if (S1L >= 0.0) {  //  ----  Region L*  ----
            EXPAND(                                                   ;  ,
                    Flux(Xt,k,j,i) = fluxL[Xt] + sl*(usL[Xt] - uL[Xt]);  ,
                    Flux(Xb,k,j,i) = fluxL[Xb] + sl*(usL[Xb] - uL[Xb]);
            )
            EXPAND(                                                       ;  ,
                    Flux(BXt,k,j,i) = fluxL[BXt] + sl*(usL[BXt] - uL[BXt]);  ,
                    Flux(BXb,k,j,i) = fluxL[BXb] + sl*(usL[BXb] - uL[BXb]);
            )
          } else if (S1R <= 0.0) { //  ----  Region R*  ----
              EXPAND(                                                   ;  ,
                      Flux(Xt,k,j,i) = fluxR[Xt] + sr*(usR[Xt] - uR[Xt]);  ,
                      Flux(Xb,k,j,i) = fluxR[Xb] + sr*(usR[Xb] - uR[Xb]);
              )
              EXPAND(                                                       ;  ,
                      Flux(BXt,k,j,i) = fluxR[BXt] + sr*(usR[BXt] - uR[BXt]);  ,
                      Flux(BXb,k,j,i) = fluxR[BXb] + sr*(usR[BXb] - uR[BXb]);
              )
          } else {
            /* ---------------------------
                  Compute U** = Uc
            --------------------------- */

            sBx = (Bx > 0.0 ? 1.0 : -1.0);

            EXPAND(                                               ,
                    usc[Xt] = 0.5*(usR[Xt] + usL[Xt]
                             + (usR[BXt] - usL[BXt])*sBx*sqrho);  ,
                    usc[Xb] = 0.5*(   usR[Xb] + usL[Xb]
                             + (usR[BXb] - usL[BXb])*sBx*sqrho);  )

            EXPAND(                                              ,
                    usc[BXt] = 0.5*(   usR[BXt] + usL[BXt]
                              + (usR[Xt] - usL[Xt])*sBx/sqrho);  ,
                    usc[BXb] = 0.5*(   usR[BXb] + usL[BXb]
                              + (usR[Xb] - usL[Xb])*sBx/sqrho);  )

            EXPAND(                                             ,
                    Flux(Xt,k,j,i) = usc[Xt]*SM - Bx*usc[BXt];  ,
                    Flux(Xb,k,j,i) = usc[Xb]*SM - Bx*usc[BXb];  )


            EXPAND(                                                  ,
                    Flux(BXt,k,j,i) = usc[BXt]*SM - Bx*usc[Xt]/rho;  ,
                    Flux(BXb,k,j,i) = usc[BXb]*SM - Bx*usc[Xb]/rho;  )
          }
        }
#endif
      }

      //6-- Compute maximum wave speed for this sweep
      cMax(k,j,i) = cmax;

      // 7-- Store the flux in the emf components
      D_EXPAND( Et(k,j,i) = st*Flux(BXt,k,j,i);  ,
                                                 ,
                Eb(k,j,i) = sb*Flux(BXb,k,j,i);  )

#if EMF_AVERAGE == UCT_CONTACT
      int s = 0;
      if (Flux(RHO,k,j,i) >  eps_UCT_CONTACT) s =  1;
      if (Flux(RHO,k,j,i) < -eps_UCT_CONTACT) s = -1;

      SV(k,j,i) = s;

#elif EMF_AVERAGE == UCT_HLL
      SL(k,j,i) = std::fmax(ZERO_F, -sl);
      SR(k,j,i) = std::fmax(ZERO_F,  sr);
#endif
  });

  idfx::popRegion();
}

#endif // HYDRO_MHDSOLVERS_HLLDMHD_HPP_
