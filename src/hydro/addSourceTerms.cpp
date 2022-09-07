// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "hydro.hpp"
#include "dataBlock.hpp"
#include "fargo.hpp"

// Add source terms
void Hydro::AddSourceTerms(real t, real dt) {
  idfx::pushRegion("Hydro::AddSourceTerms");

  IdefixArray4D<real> Uc = this->Uc;
  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  IdefixArray2D<real> fargoVelocity = data->fargo.meanVelocity;
#ifdef ISOTHERMAL
  IdefixArray3D<real> csIsoArr = this->isoSoundSpeedArray;
#endif
#if GEOMETRY == SPHERICAL
  IdefixArray1D<real> sinx2  = data->sinx2;
  IdefixArray1D<real> tanx2  = data->tanx2;
  IdefixArray1D<real> rt = data->rt;
#endif


#ifdef ISOTHERMAL
  [[maybe_unused]] real csIso = this->isoSoundSpeed;
  [[maybe_unused]] HydroModuleStatus haveIsoCs = this->haveIsoSoundSpeed;
#endif

  bool haveRotation = this->haveRotation;
  real OmegaZ = this->OmegaZ;

  // Fargo
  bool haveFargo  = data->haveFargo;

  // shearing box (only with fargo&cartesian)
  [[maybe_unused]] real sbS = this->sbS;

  if(haveUserSourceTerm) userSourceTerm(*data, t, dt);

  idefix_for("AddSourceTerms",
             data->beg[KDIR],data->end[KDIR],
             data->beg[JDIR],data->end[JDIR],
             data->beg[IDIR],data->end[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      #if GEOMETRY == CARTESIAN
        // Manually add Coriolis force in cartesian geometry. Otherwise
        // Coriolis is treated as a modification to the fluxes
        if(haveRotation) {
          Uc(MX1,k,j,i) +=   TWO_F * dt * Vc(RHO,k,j,i) * OmegaZ * Vc(VX2,k,j,i);
          Uc(MX2,k,j,i) += - TWO_F * dt * Vc(RHO,k,j,i) * OmegaZ * Vc(VX1,k,j,i);
        }
        if(haveFargo) {
          Uc(MX1,k,j,i) +=   TWO_F * dt * Vc(RHO,k,j,i) * OmegaZ * sbS * x1(i);
        }
      #endif
      // fetch fargo velocity when required
      [[maybe_unused]] real fargoV = ZERO_F;
      if(haveFargo) {
        // No source term when CARTESIAN+Fargo
        #if GEOMETRY == POLAR
          fargoV = fargoVelocity(k,i);
        #elif GEOMETRY == SPHERICAL
          fargoV = fargoVelocity(j,i);
        #endif
      }
#if GEOMETRY == CYLINDRICAL
  #if COMPONENTS == 3
      real vphi,Sm;
      vphi = Vc(iVPHI,k,j,i);
      if(haveRotation) vphi += OmegaZ*x1(i);
      Sm = Vc(RHO,k,j,i) * vphi*vphi; // Centrifugal
      // Presure (because pressure is included in the flux, additional source terms arise)
    #ifdef ISOTHERMAL
      real c2Iso;
      if(haveIsoCs == UserDefFunction) {
        c2Iso = csIsoArr(k,j,i);
        c2Iso = c2Iso*c2Iso;
      } else {
        c2Iso = csIso*csIso;
      }
      Sm += Vc(RHO,k,j,i)*c2Iso;
    #else
      Sm += Vc(PRS,k,j,i);
    #endif // ISOTHERMAL
    #if MHD==YES
      Sm -=  Vc(iBPHI,k,j,i)*Vc(iBPHI,k,j,i); // Hoop stress
      // Magnetic pressure
      Sm += HALF_F*(EXPAND( Vc(BX1,k,j,i)*Vc(BX1,k,j,i)   ,
                            +Vc(BX2,k,j,i)*Vc(BX2,k,j,i)  ,
                            +Vc(BX3,k,j,i)*Vc(BX3,k,j,i)  ));
    #endif // MHD
      Uc(MX1,k,j,i) += dt * Sm / x1(i);
  #endif // COMPONENTS

#elif GEOMETRY == POLAR
      real vphi,Sm;
      vphi = Vc(iVPHI,k,j,i) + fargoV;
      if(haveRotation) vphi += OmegaZ*x1(i);
      Sm = Vc(RHO,k,j,i) * vphi*vphi;     // Centrifugal
      // Pressure (because we're including pressure in the flux,
      // we need that to get the radial pressure gradient)
  #ifdef ISOTHERMAL
      real c2Iso;
      if(haveIsoCs == UserDefFunction) {
        c2Iso = csIsoArr(k,j,i);
        c2Iso = c2Iso*c2Iso;
      } else {
        c2Iso = csIso*csIso;
      }
      Sm += Vc(RHO,k,j,i)*c2Iso;
  #else
      Sm += Vc(PRS,k,j,i);
  #endif // ISOTHERMAL
  #if MHD==YES
      Sm -=  Vc(iBPHI,k,j,i)*Vc(iBPHI,k,j,i); // Hoop stress
      // Magnetic pressus
      Sm += HALF_F*(EXPAND( Vc(BX1,k,j,i)*Vc(BX1,k,j,i)   ,
                            +Vc(BX2,k,j,i)*Vc(BX2,k,j,i)  ,
                            +Vc(BX3,k,j,i)*Vc(BX3,k,j,i)  ));
  #endif // MHD
      Uc(MX1,k,j,i) += dt * Sm / x1(i);

#elif GEOMETRY == SPHERICAL
      real vphi,Sm,ct;
      vphi = SELECT(ZERO_F, ZERO_F, Vc(iVPHI,k,j,i))+fargoV;
      if(haveRotation) vphi += OmegaZ*x1(i)*FABS(sinx2(j));
      // Centrifugal
      Sm = Vc(RHO,k,j,i) * (EXPAND( ZERO_F, + Vc(VX2,k,j,i)*Vc(VX2,k,j,i), + vphi*vphi));
      // Pressure curvature
  #ifdef ISOTHERMAL
      real c2Iso;
      if(haveIsoCs == UserDefFunction) {
        c2Iso = csIsoArr(k,j,i);
        c2Iso = c2Iso*c2Iso;
      } else {
        c2Iso = csIso*csIso;
      }
      Sm += 2.0*Vc(RHO,k,j,i)*c2Iso;
  #else
      Sm += 2.0*Vc(PRS,k,j,i);
  #endif // ISOTHERMAL
  #if MHD == YES
      // Hoop stress
      Sm -= EXPAND( ZERO_F   ,
                    + Vc(iBTH,k,j,i)*Vc(iBTH,k,j,i)    ,
                    + Vc(iBPHI,k,j,i)*Vc(iBPHI,k,j,i)  );
      // 2* mag pressure curvature
      Sm += EXPAND( Vc(BX1,k,j,i)*Vc(BX1,k,j,i)   ,
                    +Vc(BX2,k,j,i)*Vc(BX2,k,j,i)  ,
                    +Vc(BX3,k,j,i)*Vc(BX3,k,j,i)  );
  #endif
      Uc(MX1,k,j,i) += dt*Sm/x1(i);
  #if COMPONENTS >= 2
      ct = 1.0/tanx2(j);
       // Centrifugal
      Sm = Vc(RHO,k,j,i) * (EXPAND( ZERO_F, - Vc(iVTH,k,j,i)*Vc(iVR,k,j,i), + ct*vphi*vphi));
      // Pressure curvature
  #ifdef ISOTHERMAL
      Sm += ct * c2Iso * Vc(RHO,k,j,i);
  #else
      Sm += ct * Vc(PRS,k,j,i);
  #endif
  #if MHD == YES
      // Hoop stress
      Sm += EXPAND( ZERO_F       ,
                    + Vc(iBTH,k,j,i)*Vc(iBR,k,j,i)        ,
                    - ct*Vc(iBPHI,k,j,i)*Vc(iBPHI,k,j,i)  );
      // Magnetic pressure
      Sm += HALF_F*ct*(EXPAND( Vc(BX1,k,j,i)*Vc(BX1,k,j,i)    ,
                               +Vc(BX2,k,j,i)*Vc(BX2,k,j,i)   ,
                               +Vc(BX3,k,j,i)*Vc(BX3,k,j,i))  );
  #endif
      Uc(MX2,k,j,i) += dt*Sm / rt(i);
  #endif // COMPONENTS
#endif
    }
  );

  idfx::popRegion();
}
