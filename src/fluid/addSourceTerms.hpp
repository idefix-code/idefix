// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/fluid/addSourceTerms.hpp
//
// Last modified : 10/2023
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2023)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_ADDSOURCETERMS_HPP_
#define FLUID_ADDSOURCETERMS_HPP_

#include "fluid.hpp"
#include "dataBlock.hpp"
#include "fargo.hpp"
#include "flux.hpp"
#include "input.hpp"

template<typename Phys>
struct Fluid_AddSourceTermsFunctor {
  /// @brief Functor for Add source terms
  /// @param hydro
  //*****************************************************************
  // Functor constructor
  //*****************************************************************
  explicit Fluid_AddSourceTermsFunctor(Fluid<Phys> *hydro, real dt) {
    Uc = hydro->Uc;
    Vc = hydro->Vc;
    x1 = hydro->data->x[IDIR];
    x2 = hydro->data->x[JDIR];

    this->dt = dt;
    #if GEOMETRY == SPHERICAL
      sinx2  = hydro->data->sinx2;
      tanx2  = hydro->data->tanx2;
      rt = hydro->data->rt;
    #endif
    if constexpr(Phys::isothermal) {
      eos = *(hydro->eos.get());
    }
    haveRotation = hydro->haveRotation;
    OmegaZ = hydro->OmegaZ;
    haveFargo = hydro->data->haveFargo;
    if(haveFargo) {
      fargoVelocity = hydro->data->fargo->meanVelocity;
    }
    // shearing box (only with fargo&cartesian)
    sbS = hydro->sbS;

    // Reduced velocity of light
    this->reduced_c =  hydro->reduced_c;
  }

  //*****************************************************************
  // Functor Variables
  //*****************************************************************
  IdefixArray4D<real> Uc;
  IdefixArray4D<real> Vc;
  IdefixArray1D<real> x1;
  IdefixArray1D<real> x2;
  IdefixArray3D<real> csIsoArr;

  real dt;
#if GEOMETRY == SPHERICAL
  IdefixArray1D<real> sinx2;
  IdefixArray1D<real> tanx2;
  IdefixArray1D<real> rt;
#endif

  EquationOfState eos;

  bool haveRotation;
  real OmegaZ;

  // Fargo
  bool haveFargo;
  IdefixArray2D<real> fargoVelocity;

  // shearing box (only with fargo&cartesian)
  real sbS;

  // Reduced velocity of light
  real reduced_c;

  //*****************************************************************
  // Functor Operator
  //*****************************************************************
  KOKKOS_INLINE_FUNCTION void operator() (const int k, const int j,  const int i) const {
    #if GEOMETRY == CARTESIAN
      // Manually add Coriolis force in cartesian geometry. Otherwise
      // Coriolis is treated as a modification to the fluxes
      if(haveRotation) {
        Uc(MX1,k,j,i) +=   TWO_F * dt * Vc(RHO,k,j,i) * OmegaZ * Vc(VX2,k,j,i);
        Uc(MX2,k,j,i) += - TWO_F * dt * Vc(RHO,k,j,i) * OmegaZ * Vc(VX1,k,j,i);
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

      real Fnorm2,inv_Fnorm2,Er2,f_param2,xi,Prad;
      if constexpr(Phys::radiation) {
        Fnorm2 = EXPAND(Vc(FR1,k,j,i)*Vc(FR1,k,j,i),
                        + Vc(FR2,k,j,i)*Vc(FR2,k,j,i),
                        + Vc(FR3,k,j,i)*Vc(FR3,k,j,i));
        inv_Fnorm2 = (Fnorm2 <= 1.e-100 ? ZERO_F : ONE_F / Fnorm2);
        Er2 = Vc(ER,k,j,i)*Vc(ER,k,j,i);
        f_param2 = (Er2 < 1.e-100 ? Fnorm2/(1.e-100) : Fnorm2/(Er2));
        f_param2 = FMIN(f_param2,ONE_F);
        xi  = 3.+4.*f_param2;
        xi /= 5.+2.*std::sqrt(4.-3.*f_param2);

        Prad = HALF_F*(1.-xi)*Vc(ER,k,j,i);
      }

#if GEOMETRY == CYLINDRICAL
  #if COMPONENTS == 3
      real vphi,Sm;
      vphi = Vc(iVPHI,k,j,i);
      if(haveRotation) vphi += OmegaZ*x1(i);
      if constexpr(Phys::radiation) {  // Radiation geometrical source terms
        Sm = 2.*Prad;
        Sm += HALF_F*(3.*xi-1.)*Vc(ER,k,j,i)*Vc(iFPHI,k,j,i)*Vc(iFPHI,k,j,i)*inv_Fnorm2);
        Sm *= reduced_c;
      } else {
        Sm = Vc(RHO,k,j,i) * vphi*vphi; // Centrifugal
      }
      // Presure (because pressure is included in the flux, additional source terms arise)
      if constexpr(Phys::isothermal) {
        real c2Iso = eos.GetWaveSpeed(k,j,i);
        c2Iso *= c2Iso;

        Sm += Vc(RHO,k,j,i)*c2Iso;
      } else if constexpr(Phys::pressure) {
        Sm += Vc(PRS,k,j,i);
      }  // Pressure
      if constexpr(Phys::mhd) {
        Sm -=  Vc(iBPHI,k,j,i)*Vc(iBPHI,k,j,i); // Hoop stress
        // Magnetic pressure
        Sm += HALF_F*(EXPAND( Vc(BX1,k,j,i)*Vc(BX1,k,j,i)   ,
                              +Vc(BX2,k,j,i)*Vc(BX2,k,j,i)  ,
                              +Vc(BX3,k,j,i)*Vc(BX3,k,j,i)  ));
      } // MHD
      Uc(MX1,k,j,i) += dt * Sm / x1(i);
  #endif // COMPONENTS

#elif GEOMETRY == POLAR
      real vphi,Sm;
      vphi = Vc(iVPHI,k,j,i) + fargoV;
      if(haveRotation) vphi += OmegaZ*x1(i);
      if constexpr(Phys::radiation) { // Radiation geometrical source terms
        Sm = 2.*Prad;
        Sm += HALF_F*(3.*xi-1.)*Vc(ER,k,j,i)*Vc(iFPHI,k,j,i)*Vc(iFPHI,k,j,i)*inv_Fnorm2;
        Sm *= reduced_c;
      } else {
        Sm = Vc(RHO,k,j,i) * vphi*vphi;  // Centrifugal
      }
      // Pressure (because we're including pressure in the flux,
      // we need that to get the radial pressure gradient)
      if constexpr(Phys::isothermal) {
        real c2Iso = eos.GetWaveSpeed(k,j,i);
        c2Iso *= c2Iso;
        Sm += Vc(RHO,k,j,i)*c2Iso;
      } else if constexpr(Phys::pressure) {
        Sm += Vc(PRS,k,j,i);
      } // Pressure
      if constexpr(Phys::mhd) {
        Sm -=  Vc(iBPHI,k,j,i)*Vc(iBPHI,k,j,i); // Hoop stress
        // Magnetic pressus
        Sm += HALF_F*(EXPAND( Vc(BX1,k,j,i)*Vc(BX1,k,j,i)   ,
                              +Vc(BX2,k,j,i)*Vc(BX2,k,j,i)  ,
                              +Vc(BX3,k,j,i)*Vc(BX3,k,j,i)  ));
      } // MHD
      Uc(MX1,k,j,i) += dt * Sm / x1(i);

#elif GEOMETRY == SPHERICAL
      real vphi,Sm;
      vphi = SELECT(ZERO_F, ZERO_F, Vc(iVPHI,k,j,i))+fargoV;
      if(haveRotation) vphi += OmegaZ*x1(i)*FABS(sinx2(j));
      if constexpr(Phys::radiation) {  // Radiation geometrical source terms
        Sm = EXPAND(2.*Prad,
                    + HALF_F*(3.*xi-1.)*Vc(ER,k,j,i)*Vc(iFTH,k,j,i)*Vc(iFTH,k,j,i)*inv_Fnorm2,
                    + HALF_F*(3.*xi-1.)*Vc(ER,k,j,i)*Vc(iFPHI,k,j,i)*Vc(iFPHI,k,j,i)*inv_Fnorm2);
        Sm *= reduced_c;
      } else {  // Centrifugal
        Sm = Vc(RHO,k,j,i) * (EXPAND( ZERO_F, + Vc(VX2,k,j,i)*Vc(VX2,k,j,i), + vphi*vphi));
      }
      // Pressure curvature
      [[maybe_unused]] real c2Iso{0};
      if constexpr(Phys::isothermal) {
        c2Iso = eos.GetWaveSpeed(k,j,i);
        c2Iso *= c2Iso;
        Sm += 2.0*Vc(RHO,k,j,i)*c2Iso;
      } else if constexpr(Phys::pressure) {
        Sm += 2.0*Vc(PRS,k,j,i);
      } // Pressure
      if constexpr(Phys::mhd) {
        // Hoop stress
        Sm -= EXPAND( ZERO_F   ,
                      + Vc(iBTH,k,j,i)*Vc(iBTH,k,j,i)    ,
                      + Vc(iBPHI,k,j,i)*Vc(iBPHI,k,j,i)  );
        // 2* mag pressure curvature
        Sm += EXPAND( Vc(BX1,k,j,i)*Vc(BX1,k,j,i)   ,
                      +Vc(BX2,k,j,i)*Vc(BX2,k,j,i)  ,
                      +Vc(BX3,k,j,i)*Vc(BX3,k,j,i)  );
      } //MHD
      Uc(MX1,k,j,i) += dt*Sm/rt(i);
  #if COMPONENTS >= 2
      real ct = 1.0/tanx2(j);
       // Centrifugal
      if constexpr(Phys::radiation) {
        Sm = EXPAND(ct*Prad,
                    - HALF_F*(3.*xi-1.)*Vc(ER,k,j,i)*Vc(iFTH,k,j,i)*Vc(iFR,k,j,i)*inv_Fnorm2,
                    + ct*HALF_F*(3.*xi-1.)*Vc(ER,k,j,i)*Vc(iFPHI,k,j,i)*Vc(iFPHI,k,j,i)*inv_Fnorm2);
        Sm *= reduced_c;
      } else {
        Sm = Vc(RHO,k,j,i) * (EXPAND( ZERO_F, - Vc(iVTH,k,j,i)*Vc(iVR,k,j,i), + ct*vphi*vphi));
      }
      // Pressure curvature
      if constexpr(Phys::isothermal) {
        Sm += ct * c2Iso * Vc(RHO,k,j,i);
      } else if constexpr(Phys::pressure) {
        Sm += ct * Vc(PRS,k,j,i);
      }
      if constexpr(Phys::mhd) {
        // Hoop stress
        Sm += EXPAND( ZERO_F       ,
                      + Vc(iBTH,k,j,i)*Vc(iBR,k,j,i)        ,
                      - ct*Vc(iBPHI,k,j,i)*Vc(iBPHI,k,j,i)  );
        // Magnetic pressure
        Sm += HALF_F*ct*(EXPAND( Vc(BX1,k,j,i)*Vc(BX1,k,j,i)    ,
                                +Vc(BX2,k,j,i)*Vc(BX2,k,j,i)   ,
                                +Vc(BX3,k,j,i)*Vc(BX3,k,j,i))  );
      } // MHD
      Uc(MX2,k,j,i) += dt*Sm / rt(i);
  #endif // COMPONENTS
#endif
    }
};


// Add source terms
template <typename Phys>
void Fluid<Phys>::AddSourceTerms(real t, real dt) {
  idfx::pushRegion("Fluid_"+prefix+"::AddSourceTerms");

  if(haveUserSourceTerm) {
    if(userSourceTerm != NULL) {
      userSourceTerm(this, t, dt);
    } else {
      // Deprecated version
      userSourceTermOld(*data, t, dt);
    }
  }

  auto func = Fluid_AddSourceTermsFunctor<Phys>(this,dt);

  idefix_for("AddSourceTerms",
             data->beg[KDIR],data->end[KDIR],
             data->beg[JDIR],data->end[JDIR],
             data->beg[IDIR],data->end[IDIR],
            func);

  idfx::popRegion();
}
#endif //FLUID_ADDSOURCETERMS_HPP_
