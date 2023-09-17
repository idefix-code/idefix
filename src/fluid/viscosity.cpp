// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// This source code is largely inspired from the viscous_flux of Pluto4.2
// ((c) P. Tzeferacos & A. Mignone)

#include <string>

#include "viscosity.hpp"
#include "dataBlock.hpp"
#include "fluid.hpp"
#include "fargo.hpp"


#define D_DX_I(q,n)  (q(n,k,j,i) - q(n,k,j,i - 1))
#define D_DY_J(q,n)  (q(n,k,j,i) - q(n,k,j - 1,i))
#define D_DZ_K(q,n)  (q(n,k,j,i) - q(n,k - 1,j,i))

#define D_DY_I(q,n)  (  0.25*(q(n,k,j + 1,i) + q(n,k,j + 1,i - 1)) \
                    - 0.25*(q(n,k,j - 1,i) + q(n,k,j - 1,i - 1)))

#define D_DZ_I(q,n)  (  0.25*(q(n,k + 1,j,i) + q(n,k + 1,j,i - 1))  \
                    - 0.25*(q(n,k - 1,j,i) + q(n,k - 1,j,i - 1)))

#define D_DX_J(q,n)  (  0.25*(q(n,k,j,i + 1) + q(n,k,j - 1,i + 1)) \
                    - 0.25*(q(n,k,j,i - 1) + q(n,k,j - 1,i - 1)))

#define D_DZ_J(q,n)  (  0.25*(q(n,k + 1,j,i) + q(n,k + 1,j - 1,i)) \
                    - 0.25*(q(n,k - 1,j,i) + q(n,k - 1,j - 1,i)))

#define D_DX_K(q,n)  (  0.25*(q(n,k,j,i + 1) + q(n,k - 1,j,i + 1)) \
                    - 0.25*(q(n,k,j,i - 1) + q(n,k - 1,j,i - 1)))

#define D_DY_K(q,n)  (  0.25*(q(n,k,j + 1,i) + q(n,k - 1,j + 1,i)) \
                    - 0.25*(q(n,k,j - 1,i) + q(n,k - 1,j - 1,i)))




// This function is technically part of the constructor,
// but since constructors cannot be Lambda-captured by cuda
// we need an extra public function for array initialisation in the constructor

void Viscosity::InitArrays() {
  // Allocate and fill arrays when needed
  #if GEOMETRY != CARTESIAN
    one_dmu = IdefixArray1D<real>("Viscosity_1dmu", data->np_tot[JDIR]);
    IdefixArray1D<real> dmu = one_dmu;
    IdefixArray1D<real> th = data->x[JDIR];
    idefix_for("ViscousInitGeometry",1,data->np_tot[JDIR],
      KOKKOS_LAMBDA(int j) {
        real scrch =  FABS((1.0-cos(th(j)))*(sin(th(j)) >= 0.0 ? 1.0:-1.0)
                     -(1.0-cos(th(j-1))) * (sin(th(j-1)) > 0.0 ? 1.0:-1.0));
        dmu(j) = 1.0/scrch;
      });
  #endif
  viscSrc = IdefixArray4D<real>("Viscosity_source", COMPONENTS, data->np_tot[KDIR],
                                                                data->np_tot[JDIR],
                                                                data->np_tot[IDIR]);
}
void Viscosity::ShowConfig() {
  if(status.status==Constant) {
    idfx::cout << "Viscosity: ENABLED with constant viscosity eta1="
                    << this->eta1 << " and eta2=" << this->eta2 << " ."<< std::endl;
  } else if (status.status==UserDefFunction) {
    idfx::cout << "Viscosity: ENABLED with user-defined viscosity function."
                   << std::endl;
    if(!viscousDiffusivityFunc) {
      IDEFIX_ERROR("No viscosity function has been enrolled");
    }
  } else {
    IDEFIX_ERROR("Unknown viscosity mode");
  }
  if(this->status.isExplicit) {
    idfx::cout << "Viscosity: uses an explicit time integration." << std::endl;
  } else if(this->status.isRKL) {
    idfx::cout << "Viscosity: uses a Runge-Kutta-Legendre time integration."
                << std::endl;
  } else {
    IDEFIX_ERROR("Unknown time integrator for viscosity.");
  }
}

void Viscosity::EnrollViscousDiffusivity(ViscousDiffusivityFunc myFunc) {
  if(this->status.status < UserDefFunction) {
    IDEFIX_WARNING("Viscous diffusivity enrollment requires Hydro/Viscosity "
                 "to be set to userdef in .ini file");
  }
  this->viscousDiffusivityFunc = myFunc;
}

// This function computes the viscous flux and stores it in Flux
// (this avoids an extra array)
// Associated source terms, present in non-cartesian geometry are also computed
// and stored in this->viscSrc for later use (in calcRhs).
void Viscosity::AddViscousFlux(int dir, const real t, const IdefixArray4D<real> &Flux) {
  idfx::pushRegion("Viscosity::AddViscousFlux");
  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> viscSrc = this->viscSrc;
  IdefixArray3D<real> dMax = this->dMax;
  IdefixArray3D<real> eta1Arr = this->eta1Arr;
  IdefixArray3D<real> eta2Arr = this->eta2Arr;
  IdefixArray1D<real> one_dmu = this->one_dmu;
  IdefixArray1D<real> x1 = this->data->x[IDIR];
  IdefixArray1D<real> x2 = this->data->x[JDIR];
  IdefixArray1D<real> x3 = this->data->x[KDIR];
  IdefixArray1D<real> x1l = this->data->xl[IDIR];
  IdefixArray1D<real> x2l = this->data->xl[JDIR];
  IdefixArray1D<real> x3l = this->data->xl[KDIR];
  IdefixArray1D<real> dx1 = this->data->dx[IDIR];
  IdefixArray1D<real> dx2 = this->data->dx[JDIR];
  IdefixArray1D<real> dx3 = this->data->dx[KDIR];

  // Fargo variables use to correct energy fluxes
  #if HAVE_ENERGY
    [[maybe_unused]] IdefixArray2D<real> fargoVelocity;
    [[maybe_unused]] Fargo::FargoType fargoType;
    bool haveFargo;
  #endif

  #if GEOMETRY == SPHERICAL
  IdefixArray1D<real> sinx2 = this->data->sinx2;
  IdefixArray1D<real> tanx2 = this->data->tanx2;
  IdefixArray1D<real> sinx2m = this->data->sinx2m;
  IdefixArray1D<real> tanx2m = this->data->tanx2m;
  #endif


  HydroModuleStatus haveViscosity = this->status.status;

  // Compute viscosity if needed
  if(haveViscosity == UserDefFunction && dir == IDIR) {
    if(viscousDiffusivityFunc) {
      viscousDiffusivityFunc(*data, t, eta1Arr, eta2Arr);
    } else {
      IDEFIX_ERROR("No user-defined viscosity function has been enrolled");
    }
  }

  #if HAVE_ENERGY
    haveFargo = this->data->haveFargo;
    if(haveFargo) {
      fargoVelocity = data->fargo->meanVelocity;
      fargoType = data->fargo->type;
    }
    [[maybe_unused]] real sbS = this->sbS;
  #endif


  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->data->beg[IDIR];
  iend = this->data->end[IDIR];
  jbeg = this->data->beg[JDIR];
  jend = this->data->end[JDIR];
  kbeg = this->data->beg[KDIR];
  kend = this->data->end[KDIR];


  real eta1Constant = this->eta1;
  real eta2Constant = this->eta2;

  ///////////////////////////////////////////
  // IDIR sweep                            //
  ///////////////////////////////////////////
  if(dir==IDIR) {
    iend++;
    idefix_for("ViscousFluxIDIR",kbeg, kend, jbeg, jend, ibeg, iend,
      KOKKOS_LAMBDA (int k, int j, int i) {
        [[maybe_unused]] real tau_xx, tau_xy, tau_xz;
        [[maybe_unused]] real tau_yy, tau_yz;
        [[maybe_unused]] real tau_zz;

        [[maybe_unused]] real dVxi, dVxj, dVxk;
        [[maybe_unused]] real dVyi, dVyj, dVyk;
        [[maybe_unused]] real dVzi, dVzj, dVzk;

        [[maybe_unused]] real divV;

        tau_xx = tau_xy = tau_xz = ZERO_F;
        tau_yy = tau_yz = ZERO_F;
        tau_zz = ZERO_F;

        dVxi = dVxj = dVxk = ZERO_F;
        dVyi = dVyj = dVyk = ZERO_F;
        dVzi = dVzj = dVzk = ZERO_F;

        real eta1, eta2;
        [[maybe_unused]] real etaC1, etaC2;

        if(haveViscosity == UserDefFunction) {
          etaC1 = eta1Arr(k,j ,i);
          eta1 = HALF_F*(eta1Arr(k,j,i-1)+eta1Arr(k,j,i));
          etaC2 = eta2Arr(k,j,i);
          eta2 = HALF_F*(eta2Arr(k,j,i-1)+eta2Arr(k,j,i));
        } else {
          etaC1 = eta1 = eta1Constant;
          etaC2 = eta2 = eta2Constant;
        }

        EXPAND(  dVxi = D_DX_I(Vc,VX1)/dx1(i); ,
                  dVyi = D_DX_I(Vc,VX2)/dx1(i); ,
                  dVzi = D_DX_I(Vc,VX3)/dx1(i); )

        #if DIMENSIONS >= 2
          EXPAND(  dVxj = D_DY_I(Vc,VX1)/dx2(j); ,
                    dVyj = D_DY_I(Vc,VX2)/dx2(j); ,
                    dVzj = D_DY_I(Vc,VX3)/dx2(j); )
          #if DIMENSIONS == 3
            dVxk = D_DZ_I(Vc,VX1)/dx3(k);
            dVyk = D_DZ_I(Vc,VX2)/dx3(k);
            dVzk = D_DZ_I(Vc,VX3)/dx3(k);
          #endif
        #endif

        #if GEOMETRY == CARTESIAN
          divV = D_EXPAND(dVxi, + dVyj, + dVzk);

          tau_xx = 2.0*eta1*dVxi + (eta2 - (2.0/3.0)*eta1)*divV;
          tau_xy = eta1*(dVxj + dVyi);
          tau_xz = eta1*(dVxk + dVzi);

        // No source term in cartesian geometry
        #endif // GEOMETRY == CARTESIAN
        #if GEOMETRY == CYLINDRICAL
          real one_dVr = x1(i)*FABS(x1(i))-x1(i-1)*FABS(x1(i-1));
          real drrVr = (Vc(VX1,k,j,i)*x1(i)- Vc(VX1,k,j,i-1)*FABS(x1(i-1)))*one_dVr;
          divV = D_EXPAND( drrVr, + dVyj, +ZERO_F);

          tau_xx = 2.0*eta1*drrVr + (eta2 - (2.0/3.0)*eta1)*divV;
          tau_xy = eta1*(dVxj + dVyi);
          SELECT(  ,
                    ,
                  tau_xz = eta1*0.5*(x1(i)+x1(i-1))/dx1(i)
                          *(Vc(VX3,k,j,i)/x1(i) - Vc(VX3,k,j,i-1)/x1(i-1)); )

          tau_xz = eta1*(dVxk + dVzi);

          // compute tau_zz at cell center
          divV = D_EXPAND(  0.5*(Vc(VX1,k,j,i+1) - Vc(VX1,k,j,i-1))/dx1(i) + Vc(VX1,k,j,i)/x1(i),
                            +0.5*(Vc(VX2,k,j+1,i) - Vc(VX2,k,j-1,i))/dx2(j)                      ,
                                                                                                  );
          tau_zz = 2.0*etaC1*Vc(VX1,k,j,i)/x1(i) + (etaC2 - (2.0/3.0)*etaC1)*divV;

          EXPAND( viscSrc(IDIR,k,j,i) =  -tau_zz/x1(i);  ,
                  viscSrc(JDIR,k,j,i) = ZERO_F;          ,
                  viscSrc(KDIR,k,j,i) = ZERO_F;          )


        #endif // GEOMETRY == CYLINDRICAL
        #if GEOMETRY == POLAR
          real vx1i = 0.5*(Vc(VX1,k,j,i-1)+Vc(VX1,k,j,i));
          #if COMPONENTS >= 2
            real vx2i = 0.5*(Vc(VX2,k,j,i-1)+Vc(VX2,k,j,i));
          #endif

          divV = D_EXPAND(vx1i/x1l(i) + dVxi, + dVyj/x1l(i), + dVzk);

          tau_xx = 2.0*eta1*dVxi + (eta2 - (2.0/3.0)*eta1)*divV;
          #if DIMENSIONS == 1
            tau_xy = eta1*dVyi;
          #else
            tau_xy = eta1*(dVxj/x1l(i) + dVyi - vx2i/x1l(i));
          #endif
          tau_xz = eta1*(dVxk + dVzi);


          // compute tau_yy at cell center
          divV = D_EXPAND(  0.5*(Vc(VX1,k,j,i+1) - Vc(VX1,k,j,i-1))/dx1(i) + Vc(VX1,k,j,i)/x1(i),
                            +0.5*(Vc(VX2,k,j+1,i) - Vc(VX2,k,j-1,i))/dx2(j)/x1(i)                ,
                            +0.5*(Vc(VX3,k+1,j,i) - Vc(VX3,k-1,j,i))/dx3(k)                      );

          #if DIMENSIONS == 1
            tau_yy = 2.0*etaC1*( Vc(VX1,k,j,i)/x1(i)) + (etaC2 - (2.0/3.0)*etaC1)*divV;
          #else
            tau_yy = 2.0*etaC1*( 0.5*(Vc(VX2,k,j+1,i) - Vc(VX2,k,j-1,i))/dx2(j)/x1(i)
                                + Vc(VX1,k,j,i)/x1(i))
                                + (etaC2 - (2.0/3.0)*etaC1)*divV;
          #endif
          EXPAND( viscSrc(IDIR,k,j,i) =  -tau_yy/x1(i);  ,
                  viscSrc(JDIR,k,j,i) = ZERO_F;          ,
                  viscSrc(KDIR,k,j,i) = ZERO_F;          )
        #endif // GEOMETRY == POLAR
        #if GEOMETRY == SPHERICAL
          real tan_1 = tanx2(j);
          real s_1 = sinx2(j);

          // Trick to ensure that the axis does not lead to Nans
          if(FABS(tan_1) < SMALL_NUMBER) {
            tan_1 = ZERO_F;
          } else {
            tan_1 = ONE_F/tan_1;
          }

          if(FABS(s_1) < SMALL_NUMBER) {
            s_1 = ZERO_F;
          } else {
            s_1 = ONE_F/s_1;
          }

          real vx1i = 0.5*(Vc(VX1,k,j,i-1)+Vc(VX1,k,j,i));
          #if COMPONENTS >= 2
            real vx2i = 0.5*(Vc(VX2,k,j,i-1)+Vc(VX2,k,j,i));
          #endif
          divV = D_EXPAND(2.0*vx1i/x1l(i) + dVxi,
                          + dVyj/x1l(i) + tan_1*vx2i/x1l(i),
                          + dVzk/x1l(i)*s_1 );

          tau_xx = 2.0*eta1*dVxi + (eta2 - (2.0/3.0)*eta1)*divV;

          tau_xy = dVxj/x1l(i) + dVyi;
          #if COMPONENTS >= 2
          tau_xy += - vx2i/x1l(i);
          #endif
          tau_xy *= eta1;

          tau_xz = dVxk*s_1/x1l(i) + dVzi;
          #if COMPONENTS == 3
            tau_xz += - 0.5*(Vc(VX3,k,j,i-1)+Vc(VX3,k,j,i))/x1l(i);
          #endif
          tau_xz *= eta1;

          // Compute tau_yy and tau_zz at cell center

          divV = D_EXPAND( 0.5*(Vc(VX1,k,j,i+1) - Vc(VX1,k,j,i-1))/dx1(i) + 2.0*Vc(VX1,k,j,i)/x1(i),
                            +0.5*(Vc(VX2,k,j+1,i) - Vc(VX2,k,j-1,i))/dx2(j)/x1(i)
                            +Vc(VX2,k,j,i)/x1(i)*tan_1                                        ,
                            +0.5*(Vc(VX3,k+1,j,i) - Vc(VX3,k-1,j,i))/dx3(k)/x1(i)*s_1 );

          tau_yy = Vc(VX1,k,j,i)/x1(i);
          #if DIMENSIONS >= 2
            tau_yy += 0.5*(Vc(VX2,k,j+1,i) - Vc(VX2,k,j-1,i))/dx2(j)/x1(i);
          #endif
          tau_yy = 2.0*etaC1*tau_yy + (etaC2-2.0/3.0*etaC1)*divV;

          tau_zz = Vc(VX1,k,j,i)/x1(i);
          #if COMPONENTS >= 2
            tau_zz += Vc(VX2,k,j,i)*tan_1/x1(i);
          #endif
          #if DIMENSIONS == 3
            tau_zz += 0.5*(Vc(VX3,k+1,j,i) - Vc(VX3,k-1,j,i))/dx3(k)/x1(i)*s_1;
          #endif
          tau_zz = 2.0*etaC1*tau_zz + (etaC2-2.0/3.0*etaC1)*divV;

          EXPAND( viscSrc(IDIR,k,j,i) =  -(tau_yy + tau_zz)/x1(i);  ,
                  viscSrc(JDIR,k,j,i) = ZERO_F;          ,
                  viscSrc(KDIR,k,j,i) = ZERO_F;          )

        #endif // GEOMETRY == SPHERICAL

        // Update flux with the stress tensor
        EXPAND( Flux(MX1, k, j, i) -= tau_xx; ,
                Flux(MX2, k, j, i) -= tau_xy; ,
                Flux(MX3, k, j, i) -= tau_xz; )

        #if HAVE_ENERGY
          EXPAND( real vx1 = 0.5*(Vc(VX1,k,j,i) + Vc(VX1,k,j,i-1)); ,
                  real vx2 = 0.5*(Vc(VX2,k,j,i) + Vc(VX2,k,j,i-1)); ,
                  real vx3 = 0.5*(Vc(VX3,k,j,i) + Vc(VX3,k,j,i-1)); )

          if(haveFargo) {
            // fetch fargo velocity when required and substract from vxj since
            // energy equation is computed for the deviations
            [[maybe_unused]] real meanV = 0.0;
            #if (GEOMETRY == POLAR || GEOMETRY == CARTESIAN) && DIMENSIONS >=2
              if(fargoType==Fargo::userdef) {
                meanV = HALF_F*(fargoVelocity(k,i-1)+fargoVelocity(k,i));
              } else if(fargoType==Fargo::shearingbox) {
                meanV = sbS * x1l(i);
              }
              vx2 = vx2 - meanV;
            #elif GEOMETRY == SPHERICAL && DIMENSIONS ==3
              meanV = fargoVelocity(j,i);
              vx3 = vx3 - meanV;
            #endif
          }

          Flux(ENG, k, j, i) -= EXPAND(   vx1*tau_xx  ,
                                        + vx2*tau_xy  ,
                                        + vx3*tau_xz);
        #endif

        real locdmax = (FMAX(eta1,eta2))/(0.5*(Vc(RHO,k,j,i)+Vc(RHO,k,j,i-1)));
        dMax(k,j,i) = FMAX(dMax(k,j,i),locdmax);
      });

  } else if(dir==JDIR) {
    ///////////////////////////////////////////
    // JDIR sweep                            //
    ///////////////////////////////////////////

    jend++;
    idefix_for("ViscousFluxJDIR",kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      [[maybe_unused]] real tau_xx, tau_xy, tau_xz;
      [[maybe_unused]] real tau_yy, tau_yz;
      [[maybe_unused]] real tau_zz;

      [[maybe_unused]] real dVxi, dVxj, dVxk;
      [[maybe_unused]] real dVyi, dVyj, dVyk;
      [[maybe_unused]] real dVzi, dVzj, dVzk;

      real divV;

      tau_xx = tau_xy = tau_xz = ZERO_F;
      tau_yy = tau_yz = ZERO_F;
      tau_zz = ZERO_F;

      dVxi = dVxj = dVxk = ZERO_F;
      dVyi = dVyj = dVyk = ZERO_F;
      dVzi = dVzj = dVzk = ZERO_F;

      real eta1, eta2;
      [[maybe_unused]] real etaC1, etaC2;

      if(haveViscosity == UserDefFunction) {
        etaC1 = eta1Arr(k,j,i);
        eta1 = HALF_F*(eta1Arr(k,j-1,i)+eta1Arr(k,j,i));
        etaC2 = eta2Arr(k,j,i);
        eta2 = HALF_F*(eta2Arr(k,j-1,i)+eta2Arr(k,j,i));
      } else {
        etaC1 = eta1 = eta1Constant;
        etaC2 = eta2 = eta2Constant;
      }

      EXPAND(  dVxi = D_DX_J(Vc,VX1)/dx1(i); ,
                dVyi = D_DX_J(Vc,VX2)/dx1(i); ,
                dVzi = D_DX_J(Vc,VX3)/dx1(i); )

      EXPAND(  dVxj = D_DY_J(Vc,VX1)/dx2(j); ,
                dVyj = D_DY_J(Vc,VX2)/dx2(j); ,
                dVzj = D_DY_J(Vc,VX3)/dx2(j); )

      #if DIMENSIONS == 3
        dVxk = D_DZ_J(Vc,VX1)/dx3(k);
        dVyk = D_DZ_J(Vc,VX2)/dx3(k);
        dVzk = D_DZ_J(Vc,VX3)/dx3(k);
      #endif

      #if GEOMETRY == CARTESIAN
        divV = D_EXPAND(dVxi, + dVyj, + dVzk);

        tau_xy = eta1*(dVxj+dVyi);
        tau_yy = 2.0*eta1*dVyj + (eta2 - (2.0/3.0)*eta1)*divV;
        tau_yz = eta1*(dVzj + dVyk);
      #endif // GEOMETRY == CARTESIAN

      #if GEOMETRY == CYLINDRICAL
        real vx1i = 0.5*(Vc(VX1,k,j-1,i)+Vc(VX1,k,j,i));

        divV = D_EXPAND(vx1i/x1(i) + dVxi, + dVyj, +0.0);

        tau_xy = eta1*(dVxj+dVyi);
        tau_yy = 2.0*eta1*dVyj + (eta2 - (2.0/3.0)*eta1)*divV;
        #if DIMENSIONS == 3
        tau_yz = eta1*(dVyk);
        #endif

        // no source term
        EXPAND( viscSrc(IDIR,k,j,i) = ZERO_F;  ,
                viscSrc(JDIR,k,j,i) = ZERO_F;  ,
                viscSrc(KDIR,k,j,i) = ZERO_F;  )
      #endif // GEOMETRY == CYLINDRICAL

      #if GEOMETRY == POLAR
        real vx1i = 0.5*(Vc(VX1,k,j-1,i)+Vc(VX1,k,j,i));
        real vx2i = 0.5*(Vc(VX2,k,j-1,i)+Vc(VX2,k,j,i));

        divV = D_EXPAND(vx1i/x1(i) + dVxi, + dVyj/x1(i), + dVzk);

        tau_xy = eta1*(dVxj/x1(i)+dVyi - vx2i/x1(i));
        tau_yy = 2.0*eta1*(dVyj/x1(i) + vx1i/x1(i)) + (eta2 - (2.0/3.0)*eta1)*divV;
        tau_yz = eta1*(dVzj/x1(i) + dVyk);

        // no source term
        EXPAND( viscSrc(IDIR,k,j,i) = ZERO_F;  ,
                viscSrc(JDIR,k,j,i) = ZERO_F;  ,
                viscSrc(KDIR,k,j,i) = ZERO_F;  )
      #endif

      #if GEOMETRY == SPHERICAL
        real tan_1 = tanx2m(j);
        real s_1 = sinx2m(j);

        // Trick to ensure that the axis does not lead to Nans
        if(FABS(tan_1) < SMALL_NUMBER) {
          tan_1 = ZERO_F;
        } else {
          tan_1 = ONE_F/tan_1;
        }

        if(FABS(s_1) < SMALL_NUMBER) {
          s_1 = ZERO_F;
        } else {
          s_1 = ONE_F/s_1;
        }


        real vx1i = 0.5*(Vc(VX1,k,j-1,i)+Vc(VX1,k,j,i));
        real vx2i = 0.5*(Vc(VX2,k,j-1,i)+Vc(VX2,k,j,i));

        divV = D_EXPAND( 2.0*vx1i/x1(i) + dVxi,
                        + (sinx2(j)*Vc(VX2,k,j,i) - FABS(sinx2(j-1))*Vc(VX2,k,j-1,i))/x1(i)
                          *one_dmu(j) ,
                        + dVzk/x1(i)*s_1 );

        // tau_xy is initially cell centered since it is involved in the source term
        tau_xy = etaC1*( 0.5*(Vc(VX1,k,j+1,i)-Vc(VX1,k,j-1,i))/x1(i)/dx2(j)
                        +0.5*(Vc(VX2,k,j,i+1)-Vc(VX2,k,j,i-1))/dx1(i)
                            - Vc(VX2,k,j,i)/x1(i));
        tau_yy = 2.0*eta1*(dVyj/x1(i) + vx1i/x1(i)) + (eta2 - (2.0/3.0)*eta1)*divV;

        tau_yz = dVzj/x1(i);
        #if COMPONENTS == 3
          tau_yz += -tan_1/x1(i)*0.5*(Vc(VX3,k,j-1,i)+Vc(VX3,k,j,i));
        #endif
        #if DIMENSIONS == 3
          tau_yz += s_1/x1(i)*dVyk;
        #endif
        tau_yz *= eta1;

        // Compute cell-centered divV & sources terms
        tan_1 = tanx2(j);
        s_1 = sinx2(j);

        // Trick to ensure that the axis does not lead to Nans
        if(FABS(tan_1) < SMALL_NUMBER) {
          tan_1 = ZERO_F;
        } else {
          tan_1 = ONE_F/tan_1;
        }

        if(FABS(s_1) < SMALL_NUMBER) {
          s_1 = ZERO_F;
        } else {
          s_1 = ONE_F/s_1;
        }

        divV = D_EXPAND( 0.5*(Vc(VX1,k,j,i+1) - Vc(VX1,k,j,i-1))/dx1(i) + 2.0*Vc(VX1,k,j,i)/x1(i),
                          +0.5*(Vc(VX2,k,j+1,i) - Vc(VX2,k,j-1,i))/dx2(j)/x1(i)
                          +Vc(VX2,k,j,i)/x1(i)*tan_1                                              ,
                          +0.5*(Vc(VX3,k+1,j,i) - Vc(VX3,k-1,j,i))/dx3(k)/x1(i)*s_1 );

        tau_zz = Vc(VX1,k,j,i)/x1(i);
        #if COMPONENTS >= 2
          tau_zz += Vc(VX2,k,j,i)/x1(i)*tan_1;
        #endif
        #if DIMENSIONS == 3
          tau_zz += 0.5*(Vc(VX3,k+1,j,i) - Vc(VX3,k-1,j,i))/dx3(k)/x1(i)*s_1;
        #endif
        tau_zz = 2*eta1*tau_zz + (eta2-2.0/3.0*eta1)*divV;

        EXPAND( viscSrc(IDIR,k,j,i) = ZERO_F;  ,
                viscSrc(JDIR,k,j,i) = (tau_xy - tau_zz*tan_1)/x1(i);  ,
                viscSrc(KDIR,k,j,i) = ZERO_F;  )


        // New tau_xy, this time face-centered
        tau_xy = eta1*(dVxj/x1(i)+dVyi - vx2i/x1(i));


      #endif
      // Update flux with the stress tensor
      EXPAND( Flux(MX1, k, j, i) -= tau_xy; ,
              Flux(MX2, k, j, i) -= tau_yy; ,
              Flux(MX3, k, j, i) -= tau_yz; )

      #if HAVE_ENERGY

      EXPAND( real vx1 = 0.5*(Vc(VX1,k,j-1,i) + Vc(VX1,k,j,i)); ,
                  real vx2 = 0.5*(Vc(VX2,k,j-1,i) + Vc(VX2,k,j,i)); ,
                  real vx3 = 0.5*(Vc(VX3,k,j-1,i) + Vc(VX3,k,j,i)); )

          if(haveFargo) {
            // fetch fargo velocity when required and substract from vxj since
            // energy equation is computed for the deviations
            [[maybe_unused]] real meanV = 0.0;
            #if (GEOMETRY == POLAR || GEOMETRY == CARTESIAN) && DIMENSIONS >=2
              if(fargoType==Fargo::userdef) {
                meanV = fargoVelocity(k,i);
              } else if(fargoType==Fargo::shearingbox) {
                meanV = sbS * x1(i);
              }
              vx2 = vx2 - meanV;
            #elif GEOMETRY == SPHERICAL && DIMENSIONS ==3
              meanV = HALF_F*(fargoVelocity(j,i)+fargoVelocity(j-1,i));
              vx3 = vx3 - meanV;
            #endif
          }

          Flux(ENG, k, j, i) -= EXPAND(   vx1*tau_xy  ,
                                        + vx2*tau_yy  ,
                                        + vx3*tau_yz);
      #endif

      real locdmax = (FMAX(eta1,eta2))/(0.5*(Vc(RHO,k,j,i)+Vc(RHO,k,j-1,i)));
      dMax(k,j,i) = FMAX(dMax(k,j,i),locdmax);
    });

  } else if(dir==KDIR) {
    kend++;

    ///////////////////////////////////////////
    // KDIR sweep                            //
    ///////////////////////////////////////////
    idefix_for("ViscousFluxKDIR",kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      [[maybe_unused]] real tau_xx, tau_xy, tau_xz;
      [[maybe_unused]] real tau_yy, tau_yz;
      [[maybe_unused]] real tau_zz;

      [[maybe_unused]] real dVxi, dVxj, dVxk;
      [[maybe_unused]] real dVyi, dVyj, dVyk;
      [[maybe_unused]] real dVzi, dVzj, dVzk;

      real divV;

      tau_xx = tau_xy = tau_xz = ZERO_F;
      tau_yy = tau_yz = ZERO_F;
      tau_zz = ZERO_F;

      dVxi = dVxj = dVxk = ZERO_F;
      dVyi = dVyj = dVyk = ZERO_F;
      dVzi = dVzj = dVzk = ZERO_F;

      real eta1, eta2;
      [[maybe_unused]]  real etaC1, etaC2;

      if(haveViscosity == UserDefFunction) {
        etaC1 = eta1Arr(k,j,i);
        eta1 = HALF_F*(eta1Arr(k-1,j,i)+eta1Arr(k,j,i));
        etaC2 = eta2Arr(k,j,i);
        eta2 = HALF_F*(eta2Arr(k-1,j,i)+eta2Arr(k,j,i));
      } else {
        etaC1 = eta1 = eta1Constant;
        etaC2 = eta2 = eta2Constant;
      }

      dVxi = D_DX_K(Vc,VX1)/dx1(i);
      dVyi = D_DX_K(Vc,VX2)/dx1(i);
      dVzi = D_DX_K(Vc,VX3)/dx1(i);
      dVxj = D_DY_K(Vc,VX1)/dx2(j);
      dVyj = D_DY_K(Vc,VX2)/dx2(j);
      dVzj = D_DY_K(Vc,VX3)/dx2(j);
      dVxk = D_DZ_K(Vc,VX1)/dx3(k);
      dVyk = D_DZ_K(Vc,VX2)/dx3(k);
      dVzk = D_DZ_K(Vc,VX3)/dx3(k);

      #if GEOMETRY == CARTESIAN
        divV = dVxi + dVyj + dVzk;

        tau_xz = eta1*(dVxk+dVzi);
        tau_yz = eta1*(dVyk+dVzj);
        tau_zz = 2.0*eta1*dVzk + (eta2 - (2.0/3.0)*eta1)*divV;
      #endif // GEOMETRY == CARTESIAN

      #if GEOMETRY == POLAR
        real vx1i = 0.5*(Vc(VX1,k-1,j,i)+Vc(VX1,k,j,i));

        divV = vx1i/x1(i) + dVxi + dVyj/x1(i) + dVzk;

        tau_xz = eta1*(dVxk+dVzi);
        tau_yz = eta1*(dVyk+dVzj);
        tau_zz = 2.0*eta1*dVzk + (eta2 - (2.0/3.0)*eta1)*divV;

        viscSrc(IDIR,k,j,i) = ZERO_F;
        viscSrc(JDIR,k,j,i) = ZERO_F;
        viscSrc(KDIR,k,j,i) = ZERO_F;

      #endif // GEOMETRY == POLAR

      #if GEOMETRY == SPHERICAL
        real tan_1 = tanx2(j);
        real s_1 = sinx2(j);

        // Trick to ensure that the axis does not lead to Nans
        if(FABS(tan_1) < SMALL_NUMBER) {
          tan_1 = ZERO_F;
        } else {
          tan_1 = ONE_F/tan_1;
        }

        if(FABS(s_1) < SMALL_NUMBER) {
          s_1 = ZERO_F;
        } else {
          s_1 = ONE_F/s_1;
        }


        real vx1i = 0.5*(Vc(VX1,k-1,j,i)+Vc(VX1,k,j,i));
        real vx2i = 0.5*(Vc(VX2,k-1,j,i)+Vc(VX2,k,j,i));
        real vx3i = 0.5*(Vc(VX3,k-1,j,i)+Vc(VX3,k,j,i));

        divV = 2.0*vx1i/x1(i) + dVxi + dVyj/x1(i) + tan_1*vx2i/x1(i) + dVzk/x1(i)*s_1;

        tau_xz = eta1*(dVxk*s_1/x1(i) + dVzi - vx3i/x1(i));
        tau_yz = eta1*(dVyk*s_1/x1(i) + dVzj/x1(i) - vx3i*tan_1/x1(i));
        tau_zz = 2.0*eta1*( dVzk*s_1/x1(i) + vx1i/x1(i) + vx2i*tan_1/x1(i))
                  + (eta2 - (2.0/3.0)*eta1)*divV;

        viscSrc(IDIR,k,j,i) = ZERO_F;
        viscSrc(JDIR,k,j,i) = ZERO_F;
        viscSrc(KDIR,k,j,i) = ZERO_F;

      #endif //GEOMETRY == SPHERICAL

      // Update flux with the stress tensor
      EXPAND( Flux(MX1, k, j, i) -= tau_xz; ,
              Flux(MX2, k, j, i) -= tau_yz; ,
              Flux(MX3, k, j, i) -= tau_zz; )

      #if HAVE_ENERGY

        EXPAND( real vx1 = 0.5*(Vc(VX1,k-1,j,i) + Vc(VX1,k,j,i)); ,
                real vx2 = 0.5*(Vc(VX2,k-1,j,i) + Vc(VX2,k,j,i)); ,
                real vx3 = 0.5*(Vc(VX3,k-1,j,i) + Vc(VX3,k,j,i)); )

          if(haveFargo) {
            // fetch fargo velocity when required and substract from vxj since
            // energy equation is computed for the deviations
            [[maybe_unused]] real meanV = 0.0;
            #if (GEOMETRY == POLAR || GEOMETRY == CARTESIAN) && DIMENSIONS >=2
              if(fargoType==Fargo::userdef) {
                meanV = HALF_F*(fargoVelocity(k-1,i)+fargoVelocity(k,i));
              } else if(fargoType==Fargo::shearingbox) {
                meanV = sbS * x1l(i);
              }
              vx2 = vx2 - meanV;
            #elif GEOMETRY == SPHERICAL && DIMENSIONS ==3
              meanV = fargoVelocity(j,i);
              vx3 = vx3 - meanV;
            #endif
          }

          Flux(ENG, k, j, i) -= EXPAND(   vx1*tau_xz  ,
                                        + vx2*tau_yz  ,
                                        + vx3*tau_zz);
      #endif

      real locdmax = (FMAX(eta1,eta2))/(0.5*(Vc(RHO,k,j,i)+Vc(RHO,k-1,j,i)));
      dMax(k,j,i) = FMAX(dMax(k,j,i),locdmax);;
    });
  }

  idfx::popRegion();
}
