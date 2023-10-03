// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_RIEMANNSOLVER_MHDSOLVERS_TVDLFMHD_HPP_
#define FLUID_RIEMANNSOLVER_MHDSOLVERS_TVDLFMHD_HPP_

#include "../idefix.hpp"
#include "extrapolateToFaces.hpp"
#include "flux.hpp"
#include "convertConsToPrim.hpp"
#include "storeFlux.hpp"
#include "constrainedTransport.hpp"

// Compute Riemann fluxes from states using TVDLF solver
template <typename Phys>
template<const int DIR>
void RiemannSolver<Phys>::TvdlfMHD(IdefixArray4D<real> &Flux) {
  idfx::pushRegion("RiemannSolver::TVDLF_MHD");

  using EMF = ConstrainedTransport<Phys>;

  constexpr int ioffset = (DIR==IDIR) ? 1 : 0;
  constexpr int joffset = (DIR==JDIR) ? 1 : 0;
  constexpr int koffset = (DIR==KDIR) ? 1 : 0;

  // extension in perp to the direction of integration, as required by CT.
  constexpr int iextend = (DIR==IDIR) ? 0 : 1;
  #if DIMENSIONS > 1
    constexpr int jextend = (DIR==JDIR) ? 0 : 1;
  #else
    constexpr int jextend = 0;
  #endif
  #if DIMENSIONS > 2
    constexpr int kextend = (DIR==KDIR) ? 0 : 1;
  #else
    constexpr int kextend = 0;
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
  // Define normal, tangent and bi-tanget indices
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
      const int Xn = DIR+MX1;
      EXPAND( const int BXn = DIR+BX1;                    ,
              const int BXt = (DIR == IDIR ? BX2 : BX1);  ,
              const int BXb = (DIR == KDIR ? BX2 : BX3);   )
      // Primitive variables
      real vL[Phys::nvar];
      real vR[Phys::nvar];
      real v[Phys::nvar];

      real uL[Phys::nvar];
      real uR[Phys::nvar];

      real fluxL[Phys::nvar];
      real fluxR[Phys::nvar];

      // Load primitive variables
      extrapol.ExtrapolatePrimVar(i, j, k, vL, vR);
      vL[BXn] = Vs(DIR,k,j,i);
      vR[BXn] = vL[BXn];
#pragma unroll
      for(int nv = 0 ; nv < Phys::nvar; nv++) {
        v[nv] = HALF_F*(vL[nv] + vR[nv]);
      }

      // Get the wave speed
      // Signal speeds
      real cRL, cmax, c2Iso;
      real gpr, Bt2, B2;

      // Init c2Isothermal (used only when isothermal approx is set)
      c2Iso = ZERO_F;

#if HAVE_ENERGY
      real gamma = eos.GetGamma(v[PRS],v[RHO]);
      gpr=gamma*v[PRS];
#else
      c2Iso = HALF_F*(eos.GetWaveSpeed(k,j,i)
                    +eos.GetWaveSpeed(k-koffset,j-joffset,i-ioffset));
      c2Iso *= c2Iso;

      gpr = c2Iso*v[RHO];
#endif
      Bt2=EXPAND( ZERO_F           ,
                  + v[BXt]*v[BXt]  ,
                  + v[BXb]*v[BXb]  );

      B2=Bt2 + v[BXn]*v[BXn];

      cRL = gpr - B2;
      cRL = cRL + B2 + std::sqrt(cRL*cRL + FOUR_F*gpr*Bt2);
      cRL = std::sqrt(HALF_F * cRL/v[RHO]);

      cmax = std::fmax(std::fabs(v[Xn]+cRL),FABS(v[Xn]-cRL));

      real sl, sr;
      sl = -cmax;
      sr = cmax;


      // 2-- Compute the conservative variables
      K_PrimToCons<Phys>(uL, vL, &eos);
      K_PrimToCons<Phys>(uR, vR, &eos);

      // 3-- Compute the left and right fluxes
      K_Flux<Phys,DIR>(fluxL, vL, uL, c2Iso);
      K_Flux<Phys,DIR>(fluxR, vR, uR, c2Iso);


      // 5-- Compute the flux from the left and right states
#pragma unroll
      for(int nv = 0 ; nv < Phys::nvar; nv++) {
        Flux(nv,k,j,i) = HALF_F*(fluxL[nv] + fluxR[nv] + cmax*(uL[nv] - uR[nv]));
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
      }
      /*else if (emfAverage==EMF::uct_hlld) {
        // We do not have the Alfven speed in the HLL solver
        K_StoreHLLD<DIR>(i,j,k,st,sb,c2Iso,sl,sr,vL,vR,uL,uR,Et,Eb,aL,aR,dL,dR);
      } */
  });

  idfx::popRegion();
}

#endif // FLUID_RIEMANNSOLVER_MHDSOLVERS_TVDLFMHD_HPP_
