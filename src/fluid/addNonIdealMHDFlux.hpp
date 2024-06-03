// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_ADDNONIDEALMHDFLUX_HPP_
#define FLUID_ADDNONIDEALMHDFLUX_HPP_

#include "fluid.hpp"
#include "dataBlock.hpp"

// Compute parabolic fluxes
template<typename Phys>
template <int dir>
void Fluid<Phys>::AddNonIdealMHDFlux(const real t) {
  idfx::pushRegion("Fluid::addNonIdealMHDFlux");

  if constexpr(Phys::mhd) {
    int ioffset,joffset,koffset;

    IdefixArray4D<real> Flux = this->FluxRiemann;
    IdefixArray4D<real> Vc   = this->Vc;
    IdefixArray4D<real> Vs   = this->Vs;
    IdefixArray3D<real> dMax = this->dMax;
    IdefixArray4D<real> J    = this->J;
    IdefixArray3D<real> etaArr = this->etaOhmic;
    IdefixArray3D<real> xAmbiArr = this->xAmbipolar;

    // these two are required to ensure that the type is captured by KOKKOS_LAMBDA
    HydroModuleStatus resistivity = resistivityStatus.status;
    HydroModuleStatus ambipolar   = ambipolarStatus.status;

    bool haveResistivity{false};
    bool haveAmbipolar{false};

    if(data->rklCycle) {
      haveResistivity = resistivityStatus.isRKL;
      haveAmbipolar = ambipolarStatus.isRKL;
    } else {
      haveResistivity = resistivityStatus.isExplicit;
      haveAmbipolar = ambipolarStatus.isExplicit;
    }

    real etaConstant = this->etaO;
    real xAConstant  = this->xA;

    ioffset=joffset=koffset=0;

    switch(dir) {
      case(IDIR):
        ioffset = 1;
        break;
      case(JDIR):
        joffset = 1;
        break;
      case(KDIR):
        koffset=1;
        break;
      default:
        IDEFIX_ERROR("Wrong direction");
    }

    // Load the diffusivity array when required
    if(resistivity == UserDefFunction && dir == IDIR) {
      if(ohmicDiffusivityFunc)
        ohmicDiffusivityFunc(*data, t, etaArr);
      else
        IDEFIX_ERROR("No user-defined Ohmic diffusivity function has been enrolled");
    }

    if(ambipolar == UserDefFunction && dir == IDIR) {
      if(ambipolarDiffusivityFunc)
        ambipolarDiffusivityFunc(*data, t, xAmbiArr);
      else
        IDEFIX_ERROR("No user-defined ambipolar diffusivity function has been enrolled");
    }

    // Note the flux follows the same sign convention as the hyperbolic flux
    // HEnce signs are reversed compared to the parabolic fluxes found in Pluto 4.3
    idefix_for("CalcParabolicFlux",
              data->beg[KDIR],data->end[KDIR]+koffset,
              data->beg[JDIR],data->end[JDIR]+joffset,
              data->beg[IDIR],data->end[IDIR]+ioffset,
      KOKKOS_LAMBDA (int k, int j, int i) {
        real Jx1, Jx2, Jx3;
        real Bx1, Bx2, Bx3;
        real eta,xA;
        real locdmax = 0.0;

        int ip1;  // Offset indieces
        [[maybe_unused]] int jp1, kp1;
        ip1=i+1;
  #if DIMENSIONS >=2
        jp1 = j+1;
  #else
        jp1=j;
  #endif
  #if DIMENSIONS == 3
        kp1 = k+1;
  #else
        kp1 = k;
  #endif


        Jx1=Jx2=Jx3=ZERO_F;

        if(resistivity == Constant)
          eta = etaConstant;

        if(ambipolar == Constant)
          xA = xAConstant;

        if(dir==IDIR) {
          EXPAND(                                              ,
                  Jx3 = AVERAGE_4D_Y(J, KDIR, k, jp1, i);      ,
                  Jx2 = AVERAGE_4D_Z(J, JDIR, kp1, j, i);
                  Jx1 = AVERAGE_4D_XYZ(J, IDIR, kp1, jp1, i);  )

          EXPAND( Bx1 = Vs(BX1s,k,j,i);                             ,
                  Bx2 = HALF_F*( Vc(BX2,k,j,i-1) + Vc(BX2,k,j,i));  ,
                  Bx3 = HALF_F*( Vc(BX3,k,j,i-1) + Vc(BX3,k,j,i));  )
          if(haveResistivity) {
            if(resistivity == UserDefFunction)
              eta = AVERAGE_3D_X(etaArr,k,j,i);

            // Do not update BX2 if BX2s is defined
            #if (DIMENSIONS < 2 && COMPONENTS >= 2)
              Flux(BX2,k,j,i) += - eta * Jx3;
            #endif
            #if (DIMENSIONS < 3 && COMPONENTS == 3)
              Flux(BX3,k,j,i) +=   eta * Jx2;
            #endif

            #if HAVE_ENERGY
              Flux(ENG,k,j,i) += EXPAND( ZERO_F            ,
                                        - Bx2 * eta * Jx3  ,
                                        + Bx3 * eta * Jx2  );
            #endif

            locdmax += eta;
          }

          if(haveAmbipolar) {
            if(ambipolar == UserDefFunction)
              xA = AVERAGE_3D_X(xAmbiArr,k,j,i);

            [[maybe_unused]] real BdotB = EXPAND( Bx1*Bx1, +Bx2*Bx2, +Bx3*Bx3);

            [[maybe_unused]] real Fx2 = -xA * BdotB * Jx3;
            [[maybe_unused]] real Fx3 = xA * BdotB * Jx2;
            #if COMPONENTS == 3
              real JdotB = Jx1 * Bx1 + Jx2 * Bx2 + Jx3 * Bx3;
              Fx2 += xA * JdotB * Bx3;
              Fx3 += -xA * JdotB * Bx2;
            #endif
            #if (DIMENSIONS < 2 && COMPONENTS >= 2)
              Flux(BX2,k,j,i) += Fx2;
            #endif
            #if (DIMENSIONS < 3 && COMPONENTS == 3)
              Flux(BX3,k,j,i) += Fx3;
            #endif

            #if HAVE_ENERGY
              Flux(ENG,k,j,i) += EXPAND( ZERO_F      ,
                                        + Bx2 * Fx2  ,
                                        + Bx3 * Fx3  );
            #endif

            locdmax += xA*BdotB;
          }
        }

        if(dir==JDIR) {
          EXPAND( Jx3 = AVERAGE_4D_X(J, KDIR, k, j, ip1);   ,
                                                            ,
                  Jx1 = AVERAGE_4D_Z(J, IDIR, kp1, j, i);
                  Jx2 = AVERAGE_4D_XYZ(J, JDIR, kp1, j, ip1);  )

          EXPAND( Bx1 = HALF_F*( Vc(BX1,k,j-1,i) + Vc(BX1,k,j,i));  ,
                  Bx2 = Vs(BX2s,k,j,i);                             ,
                  Bx3 = HALF_F*( Vc(BX3,k,j-1,i) + Vc(BX3,k,j,i));  )

          if(haveResistivity) {
            if(resistivity == UserDefFunction)
              eta = AVERAGE_3D_Y(etaArr,k,j,i);

            // This term is always overwritten by CT, since this sweep is performed whenver
            // DIMENSIONS>=2
            //#if DIMENSIONS == 1
            //  Flux(BX1,k,j,i) += eta * Jx3;
            //#endif
            #if (DIMENSIONS < 3 && COMPONENTS == 3)
              Flux(BX3,k,j,i) += - eta * Jx1;
            #endif

            #if HAVE_ENERGY
              Flux(ENG,k,j,i) += EXPAND( Bx1 * eta * Jx3   ,
                                                          ,
                                        - Bx3 * eta * Jx1  );
            #endif
            locdmax += eta;
          }


          if(haveAmbipolar) {
            if(ambipolar == UserDefFunction)
              xA = AVERAGE_3D_Y(xAmbiArr,k,j,i);

            [[maybe_unused]] real BdotB = EXPAND( Bx1*Bx1, +Bx2*Bx2, +Bx3*Bx3);

            [[maybe_unused]] real Fx1 = xA * BdotB * Jx3;
            [[maybe_unused]] real Fx3 = -xA * BdotB * Jx1;
            #if COMPONENTS == 3
              real JdotB = Jx1 * Bx1 + Jx2 * Bx2 + Jx3 * Bx3;
              Fx1 += -xA * JdotB * Bx3;
              Fx3 += xA * JdotB * Bx1;
            #endif

            // This term is always overwritten by CT, since this sweep is performed whenver
            // DIMENSIONS>=2
            //#if DIMENSIONS == 1
            //  Flux(BX1,k,j,i) += Fx1;
            //#endif
            #if (DIMENSIONS < 3 && COMPONENTS == 3)
              Flux(BX3,k,j,i) += Fx3;
            #endif

            #if HAVE_ENERGY
              Flux(ENG,k,j,i) += EXPAND( + Bx1 * Fx1  ,
                                        + ZERO_F     ,
                                        + Bx3 * Fx3  );
            #endif

            locdmax += xA*BdotB;
          }
        }

        if(dir==KDIR) {
          Jx1 = AVERAGE_4D_Y(J, IDIR, k, jp1, i);
          Jx2 = AVERAGE_4D_X(J, JDIR, k, j, ip1);
          Jx3 = AVERAGE_4D_XYZ(J, KDIR, k, jp1, ip1);

          Bx1 = HALF_F*( Vc(BX1,k-1,j,i) + Vc(BX1,k,j,i));
          Bx2 = HALF_F*( Vc(BX2,k-1,j,i) + Vc(BX2,k,j,i));
          Bx3 = Vs(BX3s,k,j,i);


          if(haveResistivity) {
            if(resistivity == UserDefFunction)
              eta = AVERAGE_3D_Z(etaArr,k,j,i);

            // This ie never needed since this is overwritten by CT
            //Flux(BX1,k,j,i) += -eta * Jx2;
            //Flux(BX2,k,j,i) += eta * Jx1;

            #if HAVE_ENERGY
              Flux(ENG,k,j,i) += - Bx1 * eta * Jx2 + Bx2 * eta * Jx1;
            #endif
            locdmax += eta;
          }

          if(haveAmbipolar) {
            if(ambipolar == UserDefFunction)
              xA = AVERAGE_3D_Z(xAmbiArr,k,j,i);

            [[maybe_unused]] real BdotB = Bx1*Bx1 + Bx2*Bx2 + Bx3*Bx3;

            [[maybe_unused]] real Fx1 = -xA * BdotB * Jx2;
            [[maybe_unused]] real Fx2 = xA * BdotB * Jx1;
            #if COMPONENTS == 3
              real JdotB = Jx1 * Bx1 + Jx2 * Bx2 + Jx3 * Bx3;
              Fx1 += xA * JdotB * Bx2;
              Fx2 += -xA * JdotB * Bx1;
            #endif
            // This is never needed since this is overwritten by CT
            //Flux(BX1,k,j,i) += Fx1;
            //Flux(BX2,k,j,i) += Fx2;

            #if HAVE_ENERGY
              Flux(ENG,k,j,i) += Bx1 * Fx1  + Bx2 * Fx2;
            #endif

            locdmax += xA*BdotB;
          }
        }
        dMax(k,j,i) = FMAX(dMax(k,j,i),locdmax);
      }
    );
  }
  idfx::popRegion();
}

#endif // FLUID_ADDNONIDEALMHDFLUX_HPP_
