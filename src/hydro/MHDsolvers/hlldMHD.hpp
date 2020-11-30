// ***********************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************************

#ifndef HYDRO_MHDSOLVERS_HLLDMHD_HPP_
#define HYDRO_MHDSOLVERS_HLLDMHD_HPP_

#include "../idefix.hpp"
#include "solversMHD.hpp"

// Compute Riemann fluxes from states using HLLD solver
template<const int DIR, const int Xn, const int Xt, const int Xb,
         const int BXn, const int BXt, const int BXb>
void HlldMHD(DataBlock & data, real gamma, real C2Iso) {
  idfx::pushRegion("HLLD_MHD");
    
  int ioffset,joffset,koffset;
  int iextend, jextend,kextend;
  ioffset=joffset=koffset=0;
  // extension in perp to the direction of integration, as required by CT.
  iextend=jextend=kextend=0;

  IdefixArray4D<real> PrimL = data.PrimL;
  IdefixArray4D<real> PrimR = data.PrimR;
  IdefixArray4D<real> Flux = data.FluxRiemann;
  IdefixArray3D<real> cMax = data.cMax;

  // References to required emf components
  IdefixArray3D<real> Eb;
  IdefixArray3D<real> Et;
  
  IdefixArray3D<int> SV;

  real gamma_m1=gamma-ONE_F;

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

      Et = data.emf.ezi;
      Eb = data.emf.eyi;
      SV = data.emf.svx;
      break;
    case(JDIR):
      joffset=1;
      D_EXPAND(
                iextend = 1;
                st = +ONE_F;  ,
                              ,
                kextend = 1;
                sb = -ONE_F;  )

      Et = data.emf.ezj;
      Eb = data.emf.exj;
      SV = data.emf.svy;
      break;
    case(KDIR):
      koffset=1;
      D_EXPAND(
                iextend = 1;
                st = -ONE_F;  ,
                jextend = 1;  ,
                sb = +ONE_F;  )

      Et = data.emf.eyk;
      Eb = data.emf.exk;
      SV = data.emf.svz;
      break;
    default:
      IDEFIX_ERROR("Wrong direction");
  }

  idefix_for("CalcRiemannFlux",
             data.beg[KDIR]-kextend,data.end[KDIR]+koffset+kextend,
             data.beg[JDIR]-jextend,data.end[JDIR]+joffset+jextend,
             data.beg[IDIR]-iextend,data.end[IDIR]+ioffset+iextend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      // Primitive variables
      real vL[NVAR];
      real vR[NVAR];

      // Conservative variables
      real uL[NVAR];
      real uR[NVAR];

      // Flux (left and right)
      real fluxL[NVAR];
      real fluxR[NVAR];

      // Signal speeds
      real cL, cR, cmax;

      // 1-- Store the primitive variables on the left, right, and averaged states
      #pragma unroll
      for(int nv = 0 ; nv < NVAR; nv++) {
        vL[nv] = PrimL(nv,k,j,i);
        vR[nv] = PrimR(nv,k,j,i);
      }

      // 2-- Get the wave speed
      real gpr, b1, b2, b3, Btmag2, Bmag2;
#if HAVE_ENERGY
      gpr = gamma*vL[PRS];
#else
      gpr = C2Iso*vL[RHO];
#endif

      // -- get total field
      b1 = b2 = b3 = 0.0;
      EXPAND ( b1 = vL[BXn];  ,
               b2 = vL[BXt];  ,
               b3 = vL[BXb];  )

      Btmag2 = b2*b2 + b3*b3;
      Bmag2  = b1*b1 + Btmag2;

      cL = gpr - Bmag2;
      cL = gpr + Bmag2 + sqrt(cL*cL + 4.0*gpr*Btmag2);
      cL = sqrt(HALF_F*cL/vL[RHO]);
      
#if HAVE_ENERGY
      gpr = gamma*vR[PRS];
#else
      gpr = C2Iso*vR[RHO];
#endif

      // -- get total field
      b1 = b2 = b3 = 0.0;
      EXPAND ( b1 = vR[BXn];  ,
               b2 = vR[BXt];  ,
               b3 = vR[BXb];  )

      Btmag2 = b2*b2 + b3*b3;
      Bmag2  = b1*b1 + Btmag2;

      cR = gpr - Bmag2;
      cR = gpr + Bmag2 + sqrt(cR*cR + 4.0*gpr*Btmag2);
      cR = sqrt(HALF_F*cR/vR[RHO]);
      
      // 4.1 
      real cminL = vL[Xn] - cL;
      real cmaxL = vL[Xn] + cL;
      
      real cminR = vR[Xn] - cR;
      real cmaxR = vR[Xn] + cR;
      
      real SL = FMIN(cminL, cminR);
      real SR = FMAX(cmaxL, cmaxR);
      
      cmax  = FMAX(FABS(SL), FABS(SR));
      
      // 2-- Compute the conservative variables
      K_PrimToCons(uL, vL, gamma_m1);
      K_PrimToCons(uR, vR, gamma_m1);

      // 3-- Compute the left and right fluxes
      #pragma unroll
      for(int nv = 0 ; nv < NVAR; nv++) {
        fluxL[nv] = uL[nv];
        fluxR[nv] = uR[nv];
      }
      
      K_Flux(fluxL, vL, fluxL, C2Iso, Xn, Xt, Xb, BXn, BXt, BXb);
      K_Flux(fluxR, vR, fluxR, C2Iso, Xn, Xt, Xb, BXn, BXt, BXb);
      
      real ptR, ptL;

#if HAVE_ENERGY
      ptL  = vL[PRS] + 0.5* ( EXPAND(vL[BX1]*vL[BX1] , + vL[BX2]*vL[BX2], + vL[BX3]*vL[BX3]) );
      ptR  = vR[PRS] + 0.5* ( EXPAND(vR[BX1]*vR[BX1] , + vR[BX2]*vR[BX2], + vR[BX3]*vR[BX3]) );
#else
      ptL  = C2Iso*vL[RHO] + 0.5* (EXPAND(vL[BX1]*vL[BX1], + vL[BX2]*vL[BX2], + vL[BX3]*vL[BX3]));
      ptR  = C2Iso*vR[RHO] + 0.5* (EXPAND(vR[BX1]*vR[BX1], + vR[BX2]*vR[BX2], + vR[BX3]*vR[BX3]));
#endif

      // 5-- Compute the flux from the left and right states
      if (SL > 0) {
        #pragma unroll
        for (int nv = 0 ; nv < NFLX; nv++) {
          Flux(nv,k,j,i) = fluxL[nv];
        }
      } else if (SR < 0) {
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
        scrh = 1.0/(SR - SL);
        Bx1  = Bx = (SR*vR[BXn] - SL*vL[BXn])*scrh; 
        sBx  = (Bx > 0.0 ? 1.0 : -1.0);

        duL  = SL - vL[Xn];
        duR  = SR - vR[Xn];

        scrh = 1.0/(duR*uR[RHO] - duL*uL[RHO]);
        SM   = (duR*uR[Xn] - duL*uL[Xn] - ptR + ptL)*scrh;

        pts  = duR*uR[RHO]*ptL - duL*uL[RHO]*ptR + 
               vL[RHO]*vR[RHO]*duR*duL*(vR[Xn]- vL[Xn]);
        pts *= scrh;

        usL[RHO] = uL[RHO]*duL/(SL - SM);
        usR[RHO] = uR[RHO]*duR/(SR - SM);

        sqrL = sqrt(usL[RHO]);
        sqrR = sqrt(usR[RHO]);

        S1L = SM - fabs(Bx)/sqrL;
        S1R = SM + fabs(Bx)/sqrR;

        /* -----------------------------------------------------------------
        3d When S1L -> SL or S1R -> SR a degeneracy occurs. 
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

        if ( (S1L - SL) <  1.e-4*(SM - SL) ) revert_to_hllc = 1;
        if ( (S1R - SR) > -1.e-4*(SR - SM) ) revert_to_hllc = 1;
            
        if (revert_to_hllc) {
          scrh = 1.0/(SR - SL);
          #pragma unroll
          for(int nv = 0 ; nv < NVAR; nv++) {
            Uhll[nv]  = SR*uR[nv] - SL*uL[nv] + fluxL[nv] - fluxR[nv];
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
          scrhL = (uL[RHO]*duL*duL - Bx*Bx)/(uL[RHO]*duL*(SL - SM) - Bx*Bx);
          scrhR = (uR[RHO]*duR*duR - Bx*Bx)/(uR[RHO]*duR*(SR - SM) - Bx*Bx);
  
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
        usL[ENG] /= SL - SM;

        scrhR  = EXPAND(vR[Xn]*Bx1, + vR[Xt]*uR[BXt], + vR[Xb]*uR[BXb]);
        scrhR -= EXPAND(     SM*Bx1, +    vsR*usR[BXt], +    wsR*usR[BXb]);
        usR[ENG] = duR*uR[ENG] - ptR*vR[Xn] + pts*SM + Bx*scrhR;
        usR[ENG] /= SR - SM;

    // 3c. Compute flux when S1L > 0 or S1R < 0

        if (S1L >= 0.0) {       //  ----  Region L*
          #pragma unroll
          for(int nv = 0 ; nv < NVAR; nv++) {
            Flux(nv,k,j,i) = fluxL[nv] + SL*(usL[nv] - uL[nv]);
          }
        } else if (S1R <= 0.0) {    //  ----  Region R*
          #pragma unroll
          for(int nv = 0 ; nv < NVAR; nv++) {
            Flux(nv,k,j,i) = fluxR[nv] + SR*(usR[nv] - uR[nv]);
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
                              + SL*(usL[nv] - uL[nv]);
              }
          } else {         //  ----  Region R**
            #pragma unroll
            for(int nv = 0 ; nv < NVAR; nv++) {
              Flux(nv,k,j,i) = fluxR[nv] + S1R*(ussr[nv]  - usR[nv])
                              + SR*(usR[nv] - uR[nv]);
            }
          }
        }  // end if (S1L < 0 S1R > 0)
#else
        real usc[NVAR];
        real rho, sqrho;
        int revert_to_hll;
        
        scrh = 1.0/(SR - SL);
        duL = SL - vL[Xn];
        duR = SR - vR[Xn];

        Bx1 = Bx = (SR*vR[BXn] - SL*vL[BXn])*scrh; 

        rho                = (uR[RHO]*duR - uL[RHO]*duL)*scrh;
        Flux(RHO,k,j,i) = (SL*uR[RHO]*duR - SR*uL[RHO]*duL)*scrh;
            
        /* ---------------------------
            compute S*
        --------------------------- */

        sqrho = sqrt(rho);

        SM  = Flux(RHO,k,j,i)/rho;
        S1L = SM - fabs(Bx)/sqrho;
        S1R = SM + fabs(Bx)/sqrho;

        /* ---------------------------------------------
            Prevent degeneracies when S1L -> SL or 
            S1R -> SR. Revert to HLL if necessary.
        --------------------------------------------- */

        revert_to_hll = 0;

        if ( (S1L - SL) <  1.e-4*(SR - SL) ) revert_to_hll = 1;
        if ( (S1R - SR) > -1.e-4*(SR - SL) ) revert_to_hll = 1;

        if (revert_to_hll) {
          scrh = 1.0/(SR - SL);
          #pragma unroll
          for(int nv = 0 ; nv < NVAR; nv++) {
            Flux(nv,k,j,i) = SL*SR*(uR[nv] - uL[nv])
                            + SR*fluxL[nv] - SL*fluxR[nv];
            Flux(nv,k,j,i) *= scrh;
          }
        } else {
          Flux(Xn,k,j,i) = (SR*fluxL[Xn] - SL*fluxR[Xn] 
                          + SR*SL*(uR[Xn] - uL[Xn]))*scrh;

          Flux(BXn,k,j,i) = SR*SL*(uR[BXn] - uL[BXn])*scrh;

      /* ---------------------------
                  Compute U*  
          --------------------------- */
          
          scrhL = 1.0/((SL - S1L)*(SL - S1R));
          scrhR = 1.0/((SR - S1L)*(SR - S1R));

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
                    Flux(Xt,k,j,i) = fluxL[Xt] + SL*(usL[Xt] - uL[Xt]);  ,
                    Flux(Xb,k,j,i) = fluxL[Xb] + SL*(usL[Xb] - uL[Xb]);  
            )
            EXPAND(                                                       ;  ,
                    Flux(BXt,k,j,i) = fluxL[BXt] + SL*(usL[BXt] - uL[BXt]);  ,
                    Flux(BXb,k,j,i) = fluxL[BXb] + SL*(usL[BXb] - uL[BXb]);  
            )
          } else if (S1R <= 0.0) { //  ----  Region R*  ----
              EXPAND(                                                   ;  ,
                      Flux(Xt,k,j,i) = fluxR[Xt] + SR*(usR[Xt] - uR[Xt]);  ,
                      Flux(Xb,k,j,i) = fluxR[Xb] + SR*(usR[Xb] - uR[Xb]);  
              )
              EXPAND(                                                       ;  ,
                      Flux(BXt,k,j,i) = fluxR[BXt] + SR*(usR[BXt] - uR[BXt]);  ,
                      Flux(BXb,k,j,i) = fluxR[BXb] + SR*(usR[BXb] - uR[BXb]);  
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
#endif
  });

  idfx::popRegion();
}

#endif // HYDRO_MHDSOLVERS_HLLDMHD_HPP_
