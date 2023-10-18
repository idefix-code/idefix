// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_RIEMANNSOLVER_MHDSOLVERS_HLLDMHD_HPP_
#define FLUID_RIEMANNSOLVER_MHDSOLVERS_HLLDMHD_HPP_

#include "../idefix.hpp"
#include "extrapolateToFaces.hpp"
#include "flux.hpp"
#include "convertConsToPrim.hpp"
#include "storeFlux.hpp"
#include "constrainedTransport.hpp"



// Compute Riemann fluxes from states using HLLD solver
template <typename Phys>
template<const int DIR>
void RiemannSolver<Phys>::HlldMHD(IdefixArray4D<real> &Flux) {
  idfx::pushRegion("RiemannSolver::HLLD_MHD");

  using EMF = ConstrainedTransport<Phys>;

  constexpr int ioffset = (DIR==IDIR) ? 1 : 0;
  constexpr int joffset = (DIR==JDIR) ? 1 : 0;
  constexpr int koffset = (DIR==KDIR) ? 1 : 0;

  int perpExtension=1;
  if (hydro->emf->averaging == EMF::uct_hll
      || hydro->emf->averaging == EMF::uct_hlld) {
        // Need two cells in the perp direction for these schemes
        perpExtension= data->nghost[DIR];
  }
  // extension in perp to the direction of integration, as required by CT.
  const int iextend = (DIR==IDIR) ? 0 : perpExtension;
  #if DIMENSIONS > 1
    const int jextend = (DIR==JDIR) ? 0 : perpExtension;
  #else
    const int jextend = 0;
  #endif
  #if DIMENSIONS > 2
    const int kextend = (DIR==KDIR) ? 0 : perpExtension;
  #else
    const int kextend = 0;
  #endif

  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Vs = this->Vs;
  IdefixArray3D<real> cMax = this->cMax;

  // Required for high order interpolations
  IdefixArray1D<real> dx = this->data->dx[DIR];

  // References to required emf components
  IdefixArray3D<real> Eb;
  IdefixArray3D<real> Et;


  const typename EMF::AveragingType emfAverage = hydro->emf->averaging;

  // Required by UCT_Contact
  IdefixArray3D<real> SV;

  // Required by UCT_HLLX
  IdefixArray3D<real> aL;
  IdefixArray3D<real> aR;
  IdefixArray3D<real> dL;
  IdefixArray3D<real> dR;

  EquationOfState eos = *(hydro->eos.get());

  ExtrapolateToFaces<Phys,DIR> extrapol = *this->GetExtrapolator<DIR>();

  // st and sb will be useful only when Hall is included
  real st = ONE_F, sb = ONE_F;

  switch(DIR) {
    case(IDIR):
      D_EXPAND(
                st = -ONE_F;  ,
                  ,

                sb = +ONE_F;  )

      Et = hydro->emf->ezi;
      Eb = hydro->emf->eyi;

      SV = hydro->emf->svx;

      aL = hydro->emf->axL;
      aR = hydro->emf->axR;

      dL = hydro->emf->dxL;
      dR = hydro->emf->dxR;

      break;
#if DIMENSIONS >= 2
    case(JDIR):
      D_EXPAND(
                st = +ONE_F;  ,
                              ,

                sb = -ONE_F;  )

      Et = hydro->emf->ezj;
      Eb = hydro->emf->exj;

      SV = hydro->emf->svy;

      aL = hydro->emf->ayL;
      aR = hydro->emf->ayR;

      dL = hydro->emf->dyL;
      dR = hydro->emf->dyR;

      break;
#endif
#if DIMENSIONS == 3
    case(KDIR):

      D_EXPAND(

                st = -ONE_F;  ,
                  ,
                sb = +ONE_F;  )

      Et = hydro->emf->eyk;
      Eb = hydro->emf->exk;

      SV = hydro->emf->svz;

      aL = hydro->emf->azL;
      aR = hydro->emf->azR;

      dL = hydro->emf->dzL;
      dR = hydro->emf->dzR;
      break;
#endif
    default:
      IDEFIX_ERROR("Wrong direction");
  }

  idefix_for("CalcRiemannFlux",
             data->beg[KDIR]-kextend,data->end[KDIR]+koffset+kextend,
             data->beg[JDIR]-jextend,data->end[JDIR]+joffset+jextend,
             data->beg[IDIR]-iextend,data->end[IDIR]+ioffset+iextend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      // Init the directions (should be in the kernel for proper optimisation by the compilers)
      EXPAND( constexpr int Xn = DIR+MX1;                    ,
              constexpr int Xt = (DIR == IDIR ? MX2 : MX1);  ,
              constexpr int Xb = (DIR == KDIR ? MX2 : MX3);  )

      EXPAND( constexpr int BXn = DIR+BX1;                    ,
              constexpr int BXt = (DIR == IDIR ? BX2 : BX1);  ,
              constexpr int BXb = (DIR == KDIR ? BX2 : BX3);   )

      // Primitive variables
      real vL[Phys::nvar];
      real vR[Phys::nvar];

      extrapol.ExtrapolatePrimVar(i, j, k, vL, vR);
      vL[BXn] = Vs(DIR,k,j,i);
      vR[BXn] = vL[BXn];

      // Conservative variables
      real uL[Phys::nvar];
      real uR[Phys::nvar];

      // Flux (left and right)
      real fluxL[Phys::nvar];
      real fluxR[Phys::nvar];

      // Signal speeds
      real cL, cR, cmax, c2Iso;

      // Init c2Isothermal (used only when isothermal approx is set)
      c2Iso = ZERO_F;

      // 2-- Get the wave speed
      real gpr, b1, b2, b3, Btmag2, Bmag2;
#if HAVE_ENERGY
      real gamma = eos.GetGamma(0.5*(vL[PRS]+vR[PRS]),0.5*(vL[RHO]+vR[RHO]));
      gpr = gamma*vL[PRS];
#else
      c2Iso = HALF_F*(eos.GetWaveSpeed(k,j,i)
                    +eos.GetWaveSpeed(k-koffset,j-joffset,i-ioffset));
      c2Iso *= c2Iso;

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
      K_PrimToCons<Phys>(uL, vL, &eos);
      K_PrimToCons<Phys>(uR, vR, &eos);

      // 3-- Compute the left and right fluxes
#pragma unroll
      for(int nv = 0 ; nv < Phys::nvar; nv++) {
        fluxL[nv] = uL[nv];
        fluxR[nv] = uR[nv];
      }

      K_Flux<Phys,DIR>(fluxL, vL, fluxL, c2Iso);
      K_Flux<Phys,DIR>(fluxR, vR, fluxR, c2Iso);

      [[maybe_unused]] int revert_to_hll = 0, revert_to_hllc = 0;

#if HAVE_ENERGY
      real ptL  = vL[PRS] + HALF_F* ( EXPAND(vL[BX1]*vL[BX1]     ,
                                        + vL[BX2]*vL[BX2]   ,
                                        + vL[BX3]*vL[BX3])  );
      real ptR  = vR[PRS] + HALF_F* ( EXPAND(vR[BX1]*vR[BX1]     ,
                                        + vR[BX2]*vR[BX2]   ,
                                        + vR[BX3]*vR[BX3])  );
#endif

      // 5-- Compute the flux from the left and right states
      if (sl > 0) {
#pragma unroll
        for (int nv = 0 ; nv < Phys::nvar; nv++) {
          Flux(nv,k,j,i) = fluxL[nv];
        }
      } else if (sr < 0) {
#pragma unroll
        for (int nv = 0 ; nv < Phys::nvar; nv++) {
          Flux(nv,k,j,i) = fluxR[nv];
        }
      } else {
        real usL[Phys::nvar];
        real usR[Phys::nvar];

        real scrh, scrhL, scrhR, duL, duR, sBx, Bx, SM, S1L, S1R;

#if HAVE_ENERGY
        real Uhll[Phys::nvar];
        real pts, sqrL, sqrR;
        [[maybe_unused]] real vsL, vsR, wsL, wsR;

        // 3c. Compute U*(L), U^*(R)
        scrh = ONE_F/(sr - sl);
        Bx = (sr*vR[BXn] - sl*vL[BXn])*scrh;
        sBx  = (Bx > 0.0 ? ONE_F : -ONE_F);

        duL  = sl - vL[Xn];
        duR  = sr - vR[Xn];

        scrh = ONE_F/(duR*uR[RHO] - duL*uL[RHO]);
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

        if ( (S1L - sl) <  1.e-4*(SM - sl) ) revert_to_hllc = 1;
        if ( (S1R - sr) > -1.e-4*(sr - SM) ) revert_to_hllc = 1;

        if (revert_to_hllc) {
          scrh = ONE_F/(sr - sl);
#pragma unroll
          for(int nv = 0 ; nv < Phys::nvar; nv++) {
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

          EXPAND( usL[BXn]  = Bx;            ,
                  usL[BXt]  = uL[BXt]*scrhL;  ,
                  usL[BXb]  = uL[BXb]*scrhL;  )

          EXPAND( usR[BXn] = Bx;            ,
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

        scrhL  = EXPAND( vL[Xn]*Bx, + vL[Xt]*uL[BXt], + vL[Xb]*uL[BXb]);
        scrhL -= EXPAND( SM*Bx,     + vsL*usL[BXt],   + wsL*usL[BXb]);
        usL[ENG]  = duL*uL[ENG] - ptL*vL[Xn] + pts*SM + Bx*scrhL;
        usL[ENG] /= sl - SM;

        scrhR  = EXPAND(vR[Xn]*Bx, + vR[Xt]*uR[BXt], + vR[Xb]*uR[BXb]);
        scrhR -= EXPAND(     SM*Bx, +    vsR*usR[BXt], +    wsR*usR[BXb]);
        usR[ENG] = duR*uR[ENG] - ptR*vR[Xn] + pts*SM + Bx*scrhR;
        usR[ENG] /= sr - SM;

    // 3c. Compute flux when S1L > 0 or S1R < 0

        if (S1L >= 0.0) {       //  ----  Region L*
#pragma unroll
          for(int nv = 0 ; nv < Phys::nvar; nv++) {
            Flux(nv,k,j,i) = fluxL[nv] + sl*(usL[nv] - uL[nv]);
          }
        } else if (S1R <= 0.0) {    //  ----  Region R*
#pragma unroll
          for(int nv = 0 ; nv < Phys::nvar; nv++) {
            Flux(nv,k,j,i) = fluxR[nv] + sr*(usR[nv] - uR[nv]);
          }
        } else {   // -- This state exists only if B_x != 0
          // Compute U**
          [[maybe_unused]]real vss, wss;
          real ussl[Phys::nvar];
          real ussr[Phys::nvar];

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

          EXPAND( ussl[BXn] = ussr[BXn] = Bx;  ,

                  ussl[BXt]  = sqrL*usR[BXt] + sqrR*usL[BXt] + sqrL*sqrR*(vsR - vsL)*sBx;
                  ussl[BXt] /= sqrL + sqrR;
                  ussr[BXt]  = ussl[BXt];       ,

                  ussl[BXb]  = sqrL*usR[BXb] + sqrR*usL[BXb] + sqrL*sqrR*(wsR - wsL)*sBx;
                  ussl[BXb] /= sqrL + sqrR;
                  ussr[BXb]  = ussl[BXb];      )

          // -- Energy jump

          scrhL  = EXPAND(SM*Bx, +  vsL*usL[BXt], +  wsL*usL[BXb]);
          scrhL -= EXPAND(SM*Bx, +  vss*ussl[BXt], +  wss*ussl[BXb]);

          scrhR  = EXPAND(SM*Bx, +  vsR*usR[BXt], +  wsR*usR[BXb]);
          scrhR -= EXPAND(SM*Bx, +  vss*ussr[BXt], +  wss*ussr[BXb]);

          ussl[ENG] = usL[ENG] - sqrL*scrhL*sBx;
          ussr[ENG] = usR[ENG] + sqrR*scrhR*sBx;


          if (SM >= 0.0) { //  ----  Region L**
#pragma unroll
            for(int nv = 0 ; nv < Phys::nvar; nv++) {
              Flux(nv,k,j,i) = fluxL[nv] + S1L*(ussl[nv]  - usL[nv])
                              + sl*(usL[nv] - uL[nv]);
              }
          } else {         //  ----  Region R**
#pragma unroll
            for(int nv = 0 ; nv < Phys::nvar; nv++) {
              Flux(nv,k,j,i) = fluxR[nv] + S1R*(ussr[nv]  - usR[nv])
                              + sr*(usR[nv] - uR[nv]);
            }
          }
        }  // end if (S1L < 0 S1R > 0)
#else // No ENERGY
        real usc[Phys::nvar];
        real rho, sqrho;

        scrh = ONE_F/(sr - sl);
        duL = sl - vL[Xn];
        duR = sr - vR[Xn];

        Bx = (sr*vR[BXn] - sl*vL[BXn])*scrh;

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

        if ( (S1L - sl) <  1.e-4*(sr - sl) ) revert_to_hll = 1;
        if ( (S1R - sr) > -1.e-4*(sr - sl) ) revert_to_hll = 1;

        if (revert_to_hll) {
          scrh = ONE_F/(sr - sl);
#pragma unroll
          for(int nv = 0 ; nv < Phys::nvar; nv++) {
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

          scrhL = ONE_F/((sl - S1L)*(sl - S1R));
          scrhR = ONE_F/((sr - S1L)*(sr - S1R));

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

            sBx = (Bx > 0.0 ? ONE_F : -ONE_F);

            EXPAND(                                               ,
                    usc[Xt] = HALF_F*(usR[Xt] + usL[Xt]
                             + (usR[BXt] - usL[BXt])*sBx*sqrho);  ,
                    usc[Xb] = HALF_F*(   usR[Xb] + usL[Xb]
                             + (usR[BXb] - usL[BXb])*sBx*sqrho);  )

            EXPAND(                                              ,
                    usc[BXt] = HALF_F*(   usR[BXt] + usL[BXt]
                              + (usR[Xt] - usL[Xt])*sBx/sqrho);  ,
                    usc[BXb] = HALF_F*(   usR[BXb] + usL[BXb]
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
      if (emfAverage==EMF::arithmetic
                || emfAverage==EMF::uct0) {
        K_StoreEMF<DIR>(i,j,k,st,sb,Flux,Et,Eb);
      } else if (emfAverage==EMF::uct_contact) {
        K_StoreContact<DIR>(i,j,k,st,sb,Flux,Et,Eb,SV);
      } else if (emfAverage==EMF::uct_hll) {
        K_StoreHLL<DIR>(i,j,k,st,sb,sl,sr,vL,vR,Et,Eb,aL,aR,dL,dR);
      } else if (emfAverage==EMF::uct_hlld) {
        K_StoreHLLD<DIR>(i,j,k,st,sb,c2Iso,sl,sr,vL,vR,uL,uR,Et,Eb,aL,aR,dL,dR);
      }
  });
  idfx::popRegion();
}

#endif // FLUID_RIEMANNSOLVER_MHDSOLVERS_HLLDMHD_HPP_
