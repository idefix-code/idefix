// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// This source code is largely inspired from the viscous_flux of Pluto4.2
// ((c) P. Tzeferacos & A. Mignone)
// 
// Implementation of monotonicity-preserving viscous flux following ZuHone et al.,
// ApJ

#include <string>

#include "viscosity.hpp"
#include "dataBlock.hpp"


#define D_DX_I(q,n)  (q(n,k,j,i) - q(n,k,j,i - 1))
#define D_DY_J(q,n)  (q(n,k,j,i) - q(n,k,j - 1,i))
#define D_DZ_K(q,n)  (q(n,k,j,i) - q(n,k - 1,j,i))

#define SL_DX(q,n,iz,iy,ix)  (q(n,iz,iy,ix) - q(n,iz,iy,ix - 1))
#define SL_DY(q,n,iz,iy,ix)  (q(n,iz,iy,ix) - q(n,iz,iy - 1,ix))
#define SL_DZ(q,n,iz,iy,ix)  (q(n,iz,iy,ix) - q(n,iz - 1,iy,ix))

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



#define BX_I  Vs(BX1s,k,j,i)
#define BY_J  Vs(BX2s,k,j,i)
#define BZ_K  Vs(BX3s,k,j,i)
#define BY_I (0.25*(Vs(BX2s,k,j,i) + Vs(BX2s,k,j + 1,i) \
             + Vs(BX2s,k,j,i - 1) + Vs(BX2s,k,j + 1,i - 1)))
#define BZ_I (0.25*(Vs(BX3s,k,j,i) + Vs(BX3s,k + 1,j,i) \
             + Vs(BX3s,k,j,i - 1) + Vs(BX3s,k + 1,j,i - 1)))
#define BX_J (0.25*(Vs(BX1s,k,j,i) + Vs(BX1s,k,j,i + 1) \
             + Vs(BX1s,k,j - 1,i) + Vs(BX1s,k,j - 1,i + 1)))
#define BZ_J (0.25*(Vs(BX3s,k,j,i) + Vs(BX3s,k + 1,j,i) \
             + Vs(BX3s,k,j - 1,i) + Vs(BX3s,k + 1,j - 1,i)))
#define BX_K (0.25*(Vs(BX1s,k,j,i) + Vs(BX1s,k,j,i + 1) \
             + Vs(BX1s,k - 1,j,i) + Vs(BX1s,k - 1,j,i + 1)))
#define BY_K (0.25*(Vs(BX2s,k,j,i) + Vs(BX2s,k,j + 1,i) \
             + Vs(BX2s,k - 1,j,i) + Vs(BX2s,k - 1,j + 1,i)))

////or what should be the exact equivalent, but some reason isn't, maybe ask Geoffroy why:
//#define   BX_I 0.5*(Vc(BX1,k,j,i) + Vc(BX1,k,j,i - 1))
//#define   BY_J 0.5*(Vc(BX2,k,j,i) + Vc(BX2,k,j - 1,i))
//#define   BZ_K 0.5*(Vc(BX3,k,j,i) + Vc(BX3,k - 1,j,i))
//#define   BY_I 0.5*(Vc(BX2,k,j,i) + Vc(BX2,k,j,i - 1))
//#define   BZ_I 0.5*(Vc(BX3,k,j,i) + Vc(BX3,k,j,i - 1))
//#define   BX_J 0.5*(Vc(BX1,k,j,i) + Vc(BX1,k,j - 1,i))
//#define   BZ_J 0.5*(Vc(BX3,k,j,i) + Vc(BX3,k,j - 1,i))
//#define   BX_K 0.5*(Vc(BX1,k,j,i) + Vc(BX1,k - 1,j,i))
//#define   BY_K 0.5*(Vc(BX2,k,j,i) + Vc(BX2,k - 1,j,i))


Viscosity::Viscosity() {
  // Default constructor
}

void Viscosity::Init(Input &input, Grid &grid, Hydro *hydroin) {
  idfx::pushRegion("Viscosity::Init");
  // Save the parent hydro object
  this->hydro = hydroin;

  if(input.CheckEntry("Hydro","viscosity")>=0) {
    if(input.Get<std::string>("Hydro","viscosity",1).compare("constant") == 0) {
      this->eta1 = input.Get<real>("Hydro","viscosity",2);
      this->eta2 = input.GetOrSet<real>("Hydro","viscosity",3, 0.0);
      this->haveViscosity = Constant;
      this->haveBraginskiiViscosity = false;
    } else if(input.Get<std::string>("Hydro","viscosity",1).compare("userdef") == 0) {
      this->haveViscosity = UserDefFunction;
      this->eta1Arr = IdefixArray3D<real>("ViscosityEta1Array",hydro->data->np_tot[KDIR],
                                                               hydro->data->np_tot[JDIR],
                                                               hydro->data->np_tot[IDIR]);
      this->eta2Arr = IdefixArray3D<real>("ViscosityEta1Array",hydro->data->np_tot[KDIR],
                                                               hydro->data->np_tot[JDIR],
                                                               hydro->data->np_tot[IDIR]);
      this->haveBraginskiiViscosity = false;
    } else if(input.Get<std::string>("Hydro","viscosity",1).compare("braginskii") == 0) {
      #if MHD == NO
        IDEFIX_ERROR("MHD is not enabled."
                     "Can only enabled Braginskii viscosity when MHD is enabled.");
      #endif
      this->haveBraginskiiViscosity = true;
      this->haveSlopeLimiter = false;
      if(input.Get<std::string>("Hydro","viscosity",2).compare("constant") == 0) {
        this->haveViscosity = Constant;
        this->etaBrag = input.Get<real>("Hydro","viscosity",3);
      } else if(input.Get<std::string>("Hydro","viscosity",2).compare("userdef") == 0) {
        this->haveViscosity = UserDefFunction;
        this->etaBragArr = IdefixArray3D<real>("ViscosityEtaBragArray",hydro->data->np_tot[KDIR],
                                                                 hydro->data->np_tot[JDIR],
                                                                 hydro->data->np_tot[IDIR]);
      } else {
        IDEFIX_ERROR("Unknown Braginskii viscosity definition in idefix.ini. "
                     "Can only be constant or userdef.");
      }
    } else if(input.Get<std::string>("Hydro", "viscosity",1).compare("braginskiilim") == 0) {
      #if MHD == NO
        IDEFIX_ERROR("MHD is not enabled."
                     "Can only enabled BraginskiiLim viscosity when MHD is enabled");
      #endif
      this->haveBraginskiiViscosity = true;
      this->haveSlopeLimiter = true;
      if(input.Get<std::string>("Hydro", "viscosity",2).compare("minmod") == 0) {
        this->slopeLimiter = minmodTh;
      }
      else if(input.Get<std::string>("Hydro", "viscosity",2).compare("vanleer") == 0) {
        this->slopeLimiter = vanLeerTh;
      }
      else if(input.Get<std::string>("Hydro", "viscosity",2).compare("mclim") == 0) {
        this->slopeLimiter = monotonizedCentralTh;
      } else {
        IDEFIX_ERROR("Unknown BraginskiiLim slope limiter definition in idefix.ini. "
                     "Can only be minmod, vanleer or mclim.");
      }
      if(input.Get<std::string>("Hydro","viscosity",3).compare("constant") == 0) {
        this->haveViscosity = Constant;
        this->etaBrag = input.Get<real>("Hydro","viscosity",4);
      } else if(input.Get<std::string>("Hydro","viscosity",3).compare("userdef") == 0) {
        this->haveViscosity = UserDefFunction;
        this->etaBragArr = IdefixArray3D<real>("ViscosityEtaBragArray",hydro->data->np_tot[KDIR],
                                                                 hydro->data->np_tot[JDIR],
                                                                 hydro->data->np_tot[IDIR]);
      } else {
        IDEFIX_ERROR("Unknown Braginskii viscosity definition in idefix.ini. "
                     "Can only be constant or userdef.");
      }
    } else {
      IDEFIX_ERROR("Unknown viscosity definition in idefix.ini. "
                   "Can only be constant, userdef, braginskii or braginskiilim.");
    }
  } else {
    IDEFIX_ERROR("I cannot create a Viscosity object without viscosity defined"
                   "in the .ini file");
  }

  // Allocate and fill arrays when needed
  #if GEOMETRY != CARTESIAN
    one_dmu = IdefixArray1D<real>("Viscosity_1dmu", hydro->data->np_tot[JDIR]);
    IdefixArray1D<real> dmu = one_dmu;
    IdefixArray1D<real> th = hydro->data->x[JDIR];
    idefix_for("ViscousInitGeometry",1,hydro->data->np_tot[JDIR],
      KOKKOS_LAMBDA(int j) {
        real scrch =  FABS((1.0-cos(th(j)))*(sin(th(j)) >= 0.0 ? 1.0:-1.0)
                     -(1.0-cos(th(j-1))) * (sin(th(j-1)) > 0.0 ? 1.0:-1.0));
        dmu(j) = 1.0/scrch;
      });
  #endif
  viscSrc = IdefixArray4D<real>("Viscosity_source", COMPONENTS, hydro->data->np_tot[KDIR],
                                                                hydro->data->np_tot[JDIR],
                                                                hydro->data->np_tot[IDIR]);
  idfx::popRegion();
}

void Viscosity::ShowConfig() {
  if(haveViscosity==Constant) {
    if(haveBraginskiiViscosity) {
      idfx::cout << "Viscosity: Braginskii viscosity ENABLED with etaBrag="
            << this->etaBrag << " ."<< std::endl;
    } else {
      idfx::cout << "Viscosity: ENABLED with constant viscosity eta1="
            << this->eta1 << " and eta2=" << this->eta2 << " ."<< std::endl;
    }
  } else if (haveViscosity==UserDefFunction) {
    if(haveBraginskiiViscosity) {
      idfx::cout << "Viscosity: Braginskii viscosity ENABLED with user-defined viscosity function."
                 << std::endl;
      if(!bragViscousDiffusivityFunc) {
        IDEFIX_ERROR("No braginski viscosity function has been enrolled");
      }
    } else {
      idfx::cout << "Viscosity: ENABLED with user-defined viscosity function."
                << std::endl;
      if(!viscousDiffusivityFunc) {
        IDEFIX_ERROR("No viscosity function has been enrolled");
      }
    }
  } else {
    IDEFIX_ERROR("Unknown viscosity mode");
  }
  if(hydro->viscosityStatus.isExplicit) {
    idfx::cout << "Viscosity: uses an explicit time integration." << std::endl;
  } else if(hydro->viscosityStatus.isRKL) {
    idfx::cout << "Viscosity: uses a Runge-Kutta-Legendre time integration."
                << std::endl;
  } else {
    IDEFIX_ERROR("Unknown time integrator for viscosity.");
  }
  if(haveSlopeLimiter) {
    idfx::cout << "Viscosity: uses a slope limiter." << std::endl;
  }
}

void Viscosity::EnrollViscousDiffusivity(ViscousDiffusivityFunc myFunc) {
  if(this->haveViscosity < UserDefFunction) {
    IDEFIX_ERROR("Viscous diffusivity enrollment requires Hydro/Viscosity "
                 "to be set to userdef in .ini file");
  }
  if(haveBraginskiiViscosity) {
    IDEFIX_ERROR("Use EnrollBragViscousDiffusivity for Braginskii viscosity");
  }
  this->viscousDiffusivityFunc = myFunc;
}

void Viscosity::EnrollBragViscousDiffusivity(BragViscousDiffusivityFunc myFunc) {
  if(this->haveViscosity < UserDefFunction) {
    IDEFIX_ERROR("Viscous diffusivity enrollment requires Hydro/Viscosity "
                 "to be set to userdef in .ini file");
  }
  if(!haveBraginskiiViscosity) {
    IDEFIX_ERROR("Use EnrollThermalConductivity for isotropic viscosity");
  }
  this->bragViscousDiffusivityFunc = myFunc;
}

// This function computes the viscous flux and stores it in hydro->fluxRiemann
// (this avoids an extra array)
// Associated source terms, present in non-cartesian geometry are also computed
// and stored in this->viscSrc for later use (in calcRhs).
void Viscosity::AddViscousFlux(int dir, const real t) {
  idfx::pushRegion("Viscosity::AddViscousFlux");
  IdefixArray4D<real> Vc = this->hydro->Vc;
  IdefixArray4D<real> Vs = this->hydro->Vs;
  IdefixArray4D<real> Flux = this->hydro->FluxRiemann;
  IdefixArray4D<real> viscSrc = this->viscSrc;
  IdefixArray3D<real> dMax = this->hydro->dMax;
  IdefixArray3D<real> eta1Arr = this->eta1Arr;
  IdefixArray3D<real> eta2Arr = this->eta2Arr;
  IdefixArray3D<real> etaBragArr = this->etaBragArr;
  IdefixArray1D<real> one_dmu = this->one_dmu;
  IdefixArray1D<real> x1 = this->hydro->data->x[IDIR];
  IdefixArray1D<real> x2 = this->hydro->data->x[JDIR];
  IdefixArray1D<real> x3 = this->hydro->data->x[KDIR];
  IdefixArray1D<real> x1l = this->hydro->data->xl[IDIR];
  IdefixArray1D<real> x2l = this->hydro->data->xl[JDIR];
  IdefixArray1D<real> x3l = this->hydro->data->xl[KDIR];
  IdefixArray1D<real> dx1 = this->hydro->data->dx[IDIR];
  IdefixArray1D<real> dx2 = this->hydro->data->dx[JDIR];
  IdefixArray1D<real> dx3 = this->hydro->data->dx[KDIR];

  #if GEOMETRY == SPHERICAL
  IdefixArray1D<real> sinx2 = this->hydro->data->sinx2;
  IdefixArray1D<real> tanx2 = this->hydro->data->tanx2;
  IdefixArray1D<real> sinx2m = this->hydro->data->sinx2m;
  IdefixArray1D<real> tanx2m = this->hydro->data->tanx2m;
  #endif


  HydroModuleStatus haveViscosity = this->haveViscosity;

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->hydro->data->beg[IDIR];
  iend = this->hydro->data->end[IDIR];
  jbeg = this->hydro->data->beg[JDIR];
  jend = this->hydro->data->end[JDIR];
  kbeg = this->hydro->data->beg[KDIR];
  kend = this->hydro->data->end[KDIR];


if(!this->haveBraginskiiViscosity) {
  // Compute viscosity if needed
  if(haveViscosity == UserDefFunction && dir == IDIR) {
    if(viscousDiffusivityFunc) {
      viscousDiffusivityFunc(*this->hydro->data, t, eta1Arr, eta2Arr);
    } else {
      IDEFIX_ERROR("No user-defined viscosity function has been enrolled");
    }
  }

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
          Flux(ENG, k, j, i) -= EXPAND(   0.5*(Vc(VX1,k,j,i) + Vc(VX1,k,j,i-1))*tau_xx  ,
                + 0.5*(Vc(VX2,k,j,i) + Vc(VX2,k,j,i-1))*tau_xy  ,
                + 0.5*(Vc(VX3,k,j,i) + Vc(VX3,k,j,i-1))*tau_xz);
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

        divV = D_EXPAND(vx1i/x1(i) + dVxi, + dVyj, + 0.0);

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
        tau_xy = etaC1*( 0.5*(Vc(VX1,k,j+1,i)-Vc(VX1,k,j-1,i))/x1(i)/dx2(j) // Il manque un terme drVtheta non?
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
        Flux(ENG, k, j, i) -= EXPAND(   0.5*(Vc(VX1,k,j,i) + Vc(VX1,k,j-1,i))*tau_xy  ,
              + 0.5*(Vc(VX2,k,j,i) + Vc(VX2,k,j-1,i))*tau_yy  ,
              + 0.5*(Vc(VX3,k,j,i) + Vc(VX3,k,j-1,i))*tau_yz);
      #endif

      real locdmax = (FMAX(eta1,eta2))/(0.5*(Vc(RHO,k,j,i)+Vc(RHO,k,j-1,i)));
      dMax(k,j,i) = FMAX(dMax(k,j,i),locdmax);
    });

  } else if(dir==KDIR) {
    ///////////////////////////////////////////
    // KDIR sweep                            //
    ///////////////////////////////////////////

    kend++;

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

        viscSrc(IDIR,k,j,i) = ZERO_F; //Pourquoi pas de source terme sur ephi ?
        viscSrc(JDIR,k,j,i) = ZERO_F;
        viscSrc(KDIR,k,j,i) = ZERO_F;

      #endif //GEOMETRY == SPHERICAL

      // Update flux with the stress tensor
      EXPAND( Flux(MX1, k, j, i) -= tau_xz; ,
        Flux(MX2, k, j, i) -= tau_yz; ,
        Flux(MX3, k, j, i) -= tau_zz; )

      #if HAVE_ENERGY
        Flux(ENG, k, j, i) -= EXPAND(   0.5*(Vc(VX1,k,j,i) + Vc(VX1,k-1,j,i))*tau_xz  ,
              + 0.5*(Vc(VX2,k,j,i) + Vc(VX2,k-1,j,i))*tau_yz  ,
              + 0.5*(Vc(VX3,k,j,i) + Vc(VX3,k-1,j,i))*tau_zz);
      #endif

      real locdmax = (FMAX(eta1,eta2))/(0.5*(Vc(RHO,k,j,i)+Vc(RHO,k-1,j,i)));
      dMax(k,j,i) = FMAX(dMax(k,j,i),locdmax);
    });
    }
  } else { //end of isotropic case.
  // Braginskii Viscosity
  // Compute viscosity if needed
  if(haveViscosity == UserDefFunction && dir == IDIR) {
    if(bragViscousDiffusivityFunc) {
      bragViscousDiffusivityFunc(*this->hydro->data, t, etaBragArr);
    } else {
      IDEFIX_ERROR("No user-defined Braginskii viscosity function has been enrolled");
    }
  }

  #if GEOMETRY == CYLINDRICAL || GEOMETRY == POLAR
    IDEFIX_ERROR("Braginskii viscosity in cylindrical and polar coordinates "
                   "has not be coded so far.");
  #endif

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->hydro->data->beg[IDIR];
  iend = this->hydro->data->end[IDIR];
  jbeg = this->hydro->data->beg[JDIR];
  jend = this->hydro->data->end[JDIR];
  kbeg = this->hydro->data->beg[KDIR];
  kend = this->hydro->data->end[KDIR];

  // Determine the offset along which we do the extrapolation
  if(dir==IDIR) iend++;
  if(dir==JDIR) jend++;
  if(dir==KDIR) kend++;

  real etaBragConstant = this->etaBrag;

  idefix_for("ViscousFlux",kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      real bi, bj, bk;
      real Bmag;
      real Pnor_par;
      real bbgradV;

      real tau_xx, tau_xy, tau_xz;
      real tau_yy, tau_yz;
      real tau_zz;

      real dVxi, dVxj, dVxk;
      real dVyi, dVyj, dVyk;
      real dVzi, dVzj, dVzk;

      real divV;

      //Cell-centered values to compute source terms
      real biC, bjC, bkC;
      real BmagC;
      real Pnor_parC;
      real bbgradVC;

      real tau_xxC, tau_xyC, tau_xzC;
      real tau_yyC, tau_yzC;
      real tau_zzC;

      real dVxiC, dVxjC, dVxkC;
      real dVyiC, dVyjC, dVykC;
      real dVziC, dVzjC, dVzkC;

      real divVC;

      bi = bj = bk = ZERO_F;
      biC = bjC = bkC = ZERO_F;

      tau_xx = tau_xy = tau_xz = ZERO_F;
      tau_yy = tau_yz = ZERO_F;
      tau_zz = ZERO_F;

      tau_xxC = tau_xyC = tau_xzC = ZERO_F;
      tau_yyC = tau_yzC = ZERO_F;
      tau_zzC = ZERO_F;

      dVxi = dVxj = dVxk = ZERO_F;
      dVyi = dVyj = dVyk = ZERO_F;
      dVzi = dVzj = dVzk = ZERO_F;

      dVxiC = dVxjC = dVxkC = ZERO_F;
      dVyiC = dVyjC = dVykC = ZERO_F;
      dVziC = dVzjC = dVzkC = ZERO_F;

      real etaBrag;
      real etaBragC;

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

        //Compute values at the center of the cells
        EXPAND(  dVxiC = HALF_F*(Vc(VX1,k,j,i + 1) - Vc(VX1,k,j,i - 1))/dx1(i); ,
                 dVyiC = HALF_F*(Vc(VX2,k,j,i + 1) - Vc(VX2,k,j,i - 1))/dx1(i); ,
                 dVziC = HALF_F*(Vc(VX3,k,j,i + 1) - Vc(VX3,k,j,i - 1))/dx1(i); )

        biC = Vc(BX1,k,j,i);
        BmagC = Vc(BX1,k,j,i)*Vc(BX1,k,j,i);
        #if DIMENSIONS >= 2
          bjC = Vc(BX2,k,j,i);
          BmagC += Vc(BX2,k,j,i)*Vc(BX2,k,j,i);
          EXPAND(  dVxjC = HALF_F*(Vc(VX1,k,j + 1,i) - Vc(VX1,k,j - 1,i))/dx2(j); ,
                   dVyjC = HALF_F*(Vc(VX2,k,j + 1,i) - Vc(VX2,k,j - 1,i))/dx2(j); ,
                   dVzjC = HALF_F*(Vc(VX3,k,j + 1,i) - Vc(VX3,k,j - 1,i))/dx2(j); )
          #if DIMENSIONS == 3
            bkC = Vc(BX3,k,j,i);
            BmagC += Vc(BX3,k,j,i)*Vc(BX3,k,j,i);
            EXPAND (  dVxkC = HALF_F*(Vc(VX1,k + 1,j,i) - Vc(VX1,k - 1,j,i))/dx3(k); ,
                      dVykC = HALF_F*(Vc(VX2,k + 1,j,i) - Vc(VX2,k - 1,j,i))/dx3(k); ,
                      dVzkC = HALF_F*(Vc(VX3,k + 1,j,i) - Vc(VX3,k - 1,j,i))/dx3(k); )
          #endif
        #endif
        if(BmagC< 0.001*SMALL_NUMBER) {
           BmagC = SQRT(BmagC) + 0.000001*SMALL_NUMBER;
        } else {
           BmagC = SQRT(BmagC);
        }
        biC /= BmagC;
        bjC /= BmagC;
        bkC /= BmagC;

        bbgradVC = D_EXPAND(
    biC*biC*dVxiC                         + bjC*biC*1./x1(i)*dVxjC
                                                                          + bkC*biC*s_1/x1(i)*dVxkC,
  + biC*bjC*(dVyiC - Vc(VX2,k,j,i)/x1(i)) + bjC*bjC*(1./x1(i)*dVyjC + Vc(VX1,k,j,i)/x1(i))
                                                                          + bkC*bjC*s_1/x1(i)*dVykC,
  + biC*bkC*(dVziC - Vc(VX3,k,j,i)/x1(i)) + bjC*bkC*(1./x1(i)*dVzjC - tan_1/x1(i)*Vc(VX3,k,j,i))
                    + bkC*bkC*(s_1/x1(i)*dVzkC + Vc(VX1,k,j,i)/x1(i) + tan_1/x1(i)*Vc(VX2,k,j,i)));

        divVC = D_EXPAND(2.*Vc(VX1,k,j,i)/x1(i) + dVxiC,
                        + dVyjC/x1(i) + tan_1*Vc(VX2,k,j,i)/x1(i),
                        + dVzkC/x1(i)*s_1 );
      #endif

      ///////////////////////////////////////////
      // IDIR sweep                            //
      ///////////////////////////////////////////
      if(dir == IDIR) {
        if(haveViscosity == UserDefFunction) {
          etaBragC = etaBragArr(k,j,i);
          if(haveSlopeLimiter) {
            etaBrag = 2.*(etaBragArr(k,j,i-1)*etaBragArr(k,j,i))/(etaBragArr(k,j,i-1)+etaBragArr(k,j,i));
          } else {
          etaBrag = HALF_F*(etaBragArr(k,j,i - 1)+etaBragArr(k,j,i));
          }
        } else {
          etaBragC = etaBrag = etaBragConstant;
        }


        EXPAND(  dVxi = D_DX_I(Vc,VX1)/dx1(i); ,
                 dVyi = D_DX_I(Vc,VX2)/dx1(i); ,
                 dVzi = D_DX_I(Vc,VX3)/dx1(i); )

        bi = BX_I;
        Bmag = BX_I*BX_I;
        #if DIMENSIONS >= 2
          bj = BY_I;
          Bmag += BY_I*BY_I;
          if (haveSlopeLimiter) {
            dVxj = slopeLimiter(SL_DY(Vc,VX1,k,j + 1,i)/dx2(j+1), SL_DY(Vc,VX1,k,j + 1,i - 1)/dx2(j+1));
            dVxj = slopeLimiter(dVxj, slopeLimiter(SL_DY(Vc,VX1,k,j,i)/dx2(j), SL_DY(Vc,VX1,k,j,i - 1)/dx2(j)));
            dVyj = slopeLimiter(SL_DY(Vc,VX2,k,j + 1,i)/dx2(j+1), SL_DY(Vc,VX2,k,j + 1,i - 1)/dx2(j+1));
            dVyj = slopeLimiter(dVyj, slopeLimiter(SL_DY(Vc,VX2,k,j,i)/dx2(j), SL_DY(Vc,VX2,k,j,i - 1)/dx2(j)));
            dVzj = slopeLimiter(SL_DY(Vc,VX3,k,j + 1,i)/dx2(j+1), SL_DY(Vc,VX3,k,j + 1,i - 1)/dx2(j+1));
            dVzj = slopeLimiter(dVzj, slopeLimiter(SL_DY(Vc,VX3,k,j,i)/dx2(j), SL_DY(Vc,VX3,k,j,i - 1)/dx2(j)));
          } else {
            EXPAND(  dVxj = D_DY_I(Vc,VX1)/dx2(j); ,
                     dVyj = D_DY_I(Vc,VX2)/dx2(j); ,
                     dVzj = D_DY_I(Vc,VX3)/dx2(j); )
          }
          #if DIMENSIONS == 3
            bk = BZ_I;
            Bmag += BZ_I*BZ_I;
            if (haveSlopeLimiter) {
              dVxk = slopeLimiter(SL_DZ(Vc,VX1,k + 1,j,i)/dx3(k+1), SL_DZ(Vc,VX1,k + 1,j,i - 1)/dx3(k+1));
              dVxk = slopeLimiter(dVxk, slopeLimiter(SL_DZ(Vc,VX1,k,j,i)/dx3(k), SL_DZ(Vc,VX1,k,j,i - 1)/dx3(k)));
              dVyk = slopeLimiter(SL_DZ(Vc,VX2,k + 1,j,i)/dx3(k+1), SL_DZ(Vc,VX2,k + 1,j,i - 1)/dx3(k+1));
              dVyk = slopeLimiter(dVyk, slopeLimiter(SL_DZ(Vc,VX2,k,j,i)/dx3(k), SL_DZ(Vc,VX2,k,j,i - 1)/dx3(k)));
              dVzk = slopeLimiter(SL_DZ(Vc,VX3,k + 1,j,i)/dx3(k+1), SL_DZ(Vc,VX3,k + 1,j,i - 1)/dx3(k+1));
              dVzk = slopeLimiter(dVzk, slopeLimiter(SL_DZ(Vc,VX3,k,j,i)/dx3(k), SL_DZ(Vc,VX3,k,j,i - 1)/dx3(k)));
            } else {
              EXPAND (  dVxk = D_DZ_I(Vc,VX1)/dx3(k); ,
                        dVyk = D_DZ_I(Vc,VX2)/dx3(k); ,
                        dVzk = D_DZ_I(Vc,VX3)/dx3(k); )
            }
          #endif
        #endif
        if(Bmag< 0.001*SMALL_NUMBER) {
          Bmag = SQRT(Bmag) + 0.000001*SMALL_NUMBER;
        } else {
          Bmag = SQRT(Bmag);
        }
        bi /= Bmag;
        bj /= Bmag;
        bk /= Bmag;

        #if GEOMETRY == CARTESIAN
          bbgradV = D_EXPAND( bi*bi*dVxi + bj*bi*dVxj + bk*bi*dVxk,
                            + bi*bj*dVyi + bj*bj*dVyj + bk*bj*dVyk,
                            + bi*bk*dVzi + bj*bk*dVzj + bk*bk*dVzk);

          divV = D_EXPAND(dVxi, + dVyj, + dVzk);

          Pnor_par = 3.0*etaBrag*(bbgradV - divV/3.0);

          tau_xx = Pnor_par*(bi*bi-1./3.);
          tau_xy = Pnor_par*bi*bj;
          tau_xz = Pnor_par*bi*bk;
        // No source term in cartesian geometry
        #endif // GEOMETRY == CARTESIAN

        #if GEOMETRY == SPHERICAL
          real vx1i, vx2i, vx3i;
          vx1i = vx2i = vx3i = ZERO_F;

          vx1i = 0.5*(Vc(VX1,k,j,i-1)+Vc(VX1,k,j,i));
          #if COMPONENTS >= 2
            vx2i = 0.5*(Vc(VX2,k,j,i-1)+Vc(VX2,k,j,i));
            #if COMPONENTS == 3
              vx3i = 0.5*(Vc(VX3,k,j,i-1)+Vc(VX3,k,j,i));
            #endif
          #endif

          // NOLINT
          bbgradV = D_EXPAND(
      bi*bi*dVxi                 + bj*bi*1./x1l(i)*dVxj                 + bk*bi*s_1/x1l(i)*dVxk,
    + bi*bj*(dVyi - vx2i/x1l(i)) + bj*bj*(1./x1l(i)*dVyj + vx1i/x1l(i)) + bk*bj*s_1/x1l(i)*dVyk,
    + bi*bk*(dVzi - vx3i/x1l(i)) + bj*bk*(1./x1l(i)*dVzj - tan_1/x1l(i)*vx3i)
                                      + bk*bk*(s_1/x1l(i)*dVzk + vx1i/x1l(i) + tan_1/x1l(i)*vx2i));

          divV = D_EXPAND(2.0*vx1i/x1l(i) + dVxi,
                          + dVyj/x1l(i) + tan_1*vx2i/x1l(i),
                          + dVzk/x1l(i)*s_1 );

          Pnor_par = 3.0*etaBrag*(bbgradV - divV/3.0);

          tau_xx = Pnor_par*(bi*bi-1./3.);
          tau_xy = Pnor_par*bi*bj;
          tau_xz = Pnor_par*bi*bk;

          //cell-centered values for source terms
          Pnor_parC = 3.0*etaBragC*(bbgradVC - divVC/3.0);

          tau_yyC = Pnor_parC*(bjC*bjC - 1./3.);
          tau_zzC = Pnor_parC*(bkC*bkC - 1./3.);

          EXPAND( viscSrc(IDIR,k,j,i) =  -(tau_yyC + tau_zzC)/x1(i);  ,
                  viscSrc(JDIR,k,j,i) = ZERO_F;          ,
                  viscSrc(KDIR,k,j,i) = ZERO_F;          )
        #endif


        // Update flux with the stress tensor
        EXPAND( Flux(MX1, k, j, i) -= tau_xx; ,
                Flux(MX2, k, j, i) -= tau_xy; ,
                Flux(MX3, k, j, i) -= tau_xz; )

        #if HAVE_ENERGY
        Flux(ENG, k, j, i) -= EXPAND(   0.5*(Vc(VX1,k,j,i) + Vc(VX1,k,j,i-1))*tau_xx  ,
                                      + 0.5*(Vc(VX2,k,j,i) + Vc(VX2,k,j,i-1))*tau_xy  ,
                                      + 0.5*(Vc(VX3,k,j,i) + Vc(VX3,k,j,i-1))*tau_xz);
        #endif

        real locdmax = etaBrag/(0.5*(Vc(RHO,k,j,i)+Vc(RHO,k,j,i-1)));
        dMax(k,j,i) = FMAX(dMax(k,j,i),locdmax);
      }


      ///////////////////////////////////////////
      // JDIR sweep                            //
      ///////////////////////////////////////////
      if(dir == JDIR) {
        if(haveViscosity == UserDefFunction) {
          etaBragC = etaBragArr(k,j,i);
          if(haveSlopeLimiter) {
            etaBrag = 2.*(etaBragArr(k,j-1,i)*etaBragArr(k,j,i))/(etaBragArr(k,j-1,i)+etaBragArr(k,j,i));
          } else {
          etaBrag = HALF_F*(etaBragArr(k,j - 1,i)+etaBragArr(k,j,i));
          }
        } else {
          etaBragC = etaBrag = etaBragConstant;
        }

        if (haveSlopeLimiter) {
          dVxi = slopeLimiter(SL_DX(Vc,VX1,k,j,i + 1)/dx1(i+1), SL_DX(Vc,VX1,k,j - 1,i + 1)/dx1(i+1));
          dVxi = slopeLimiter(dVxi, slopeLimiter(SL_DX(Vc,VX1,k,j,i)/dx1(i), SL_DX(Vc,VX1,k,j - 1,i)/dx1(i)));
          dVyi = slopeLimiter(SL_DX(Vc,VX2,k,j,i + 1)/dx1(i+1), SL_DX(Vc,VX2,k,j - 1,i + 1)/dx1(i+1));
          dVyi = slopeLimiter(dVyi, slopeLimiter(SL_DX(Vc,VX2,k,j,i)/dx1(i), SL_DX(Vc,VX2,k,j - 1,i)/dx1(i)));
          dVzi = slopeLimiter(SL_DX(Vc,VX3,k,j,i + 1)/dx1(i+1), SL_DX(Vc,VX3,k,j - 1,i + 1)/dx1(i+1));
          dVzi = slopeLimiter(dVzi, slopeLimiter(SL_DX(Vc,VX3,k,j,i)/dx1(i), SL_DX(Vc,VX3,k,j - 1,i)/dx1(i)));
        } else {
        EXPAND(  dVxi = D_DX_J(Vc,VX1)/dx1(i); ,
                 dVyi = D_DX_J(Vc,VX2)/dx1(i); ,
                 dVzi = D_DX_J(Vc,VX3)/dx1(i); )
        }

        bi = BX_J;
        Bmag = BX_J*BX_J;
        #if DIMENSIONS >= 2
          bj = BY_J;
          Bmag += BY_J*BY_J;
            EXPAND(  dVxj = D_DY_J(Vc,VX1)/dx2(j); ,
                     dVyj = D_DY_J(Vc,VX2)/dx2(j); ,
                     dVzj = D_DY_J(Vc,VX3)/dx2(j); )
          #if DIMENSIONS == 3
            bk = BZ_J;
            Bmag += BZ_J*BZ_J;
            if (haveSlopeLimiter) {
              dVxk = slopeLimiter(SL_DZ(Vc,VX1,k + 1,j,i)/dx3(k+1), SL_DZ(Vc,VX1,k + 1,j - 1,i)/dx3(k+1));
              dVxk = slopeLimiter(dVxk, slopeLimiter(SL_DZ(Vc,VX1,k,j,i)/dx3(k), SL_DZ(Vc,VX1,k,j - 1,i)/dx3(k)));
              dVyk = slopeLimiter(SL_DZ(Vc,VX2,k + 1,j,i)/dx3(k+1), SL_DZ(Vc,VX2,k + 1,j - 1,i)/dx3(k+1));
              dVyk = slopeLimiter(dVyk, slopeLimiter(SL_DZ(Vc,VX2,k,j,i)/dx3(k), SL_DZ(Vc,VX2,k,j - 1,i)/dx3(k)));
              dVzk = slopeLimiter(SL_DZ(Vc,VX3,k + 1,j,i)/dx3(k+1), SL_DZ(Vc,VX3,k + 1,j - 1,i)/dx3(k+1));
              dVzk = slopeLimiter(dVzk, slopeLimiter(SL_DZ(Vc,VX3,k,j,i)/dx3(k), SL_DZ(Vc,VX3,k,j - 1,i)/dx3(k)));
            } else {
              EXPAND (  dVxk = D_DZ_J(Vc,VX1)/dx3(k); ,
                        dVyk = D_DZ_J(Vc,VX2)/dx3(k); ,
                        dVzk = D_DZ_J(Vc,VX3)/dx3(k); )
            }
          #endif
        #endif
        if(Bmag< 0.001*SMALL_NUMBER) {
          Bmag = SQRT(Bmag) + 0.000001*SMALL_NUMBER;
        } else {
          Bmag = SQRT(Bmag);
        }
        bi /= Bmag;
        bj /= Bmag;
        bk /= Bmag;

        #if GEOMETRY == CARTESIAN
          bbgradV = D_EXPAND( bi*bi*dVxi + bj*bi*dVxj + bk*bi*dVxk,
                            + bi*bj*dVyi + bj*bj*dVyj + bk*bj*dVyk,
                            + bi*bk*dVzi + bj*bk*dVzj + bk*bk*dVzk);

          divV = D_EXPAND(dVxi, + dVyj, + dVzk);

          Pnor_par = 3.0*etaBrag*(bbgradV - divV/3.0);

          tau_xy = Pnor_par*bi*bj;
          tau_yy = Pnor_par*(bj*bj - 1./3.);
          tau_yz = Pnor_par*bj*bk;
        // No source term in cartesian geometry
        #endif // GEOMETRY == CARTESIAN


        #if GEOMETRY == SPHERICAL
          tan_1 = tanx2m(j);
          s_1 = sinx2m(j);

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

          real vx1i, vx2i, vx3i;
          vx1i = vx2i = vx3i = ZERO_F;

          vx1i = 0.5*(Vc(VX1,k,j-1,i)+Vc(VX1,k,j,i));
          #if COMPONENTS >= 2
            vx2i = 0.5*(Vc(VX2,k,j-1,i)+Vc(VX2,k,j,i));
            #if COMPONENTS == 3
              vx3i = 0.5*(Vc(VX3,k,j-1,i)+Vc(VX3,k,j,i));
            #endif
          #endif

          bbgradV = D_EXPAND(
      bi*bi*dVxi + bj*bi*1./x1(i)*dVxj + bk*bi*s_1/x1(i)*dVxk,
    + bi*bj*(dVyi - vx2i/x1(i)) + bj*bj*(1./x1(i)*dVyj + vx1i/x1(i)) + bk*bj*s_1/x1(i)*dVyk,
    + bi*bk*(dVzi - vx3i/x1(i)) + bj*bk*(1./x1(i)*dVzj - tan_1/x1(i)*vx3i)
                                        + bk*bk*(s_1/x1(i)*dVzk + vx1i/x1(i) + tan_1/x1(i)*vx2i));

          divV = D_EXPAND( 2.0*vx1i/x1(i) + dVxi,
                          +(SIN(x2(j))*Vc(VX2,k,j,i) - FABS(SIN(x2(j-1)))*Vc(VX2,k,j-1,i))/x1(i)
                           *one_dmu(j) ,
                          + dVzk/x1(i)*s_1 );

          Pnor_par = 3.0*etaBrag*(bbgradV - divV/3.0);

          tau_xy = Pnor_par*bi*bj;
          tau_yy = Pnor_par*(bj*bj - 1./3.);
          tau_yz = Pnor_par*bj*bk;

          //Cell-centered values for the source terms
          tan_1 = tanx2(j);

          // Trick to ensure that the axis does not lead to Nans
          if(FABS(tan_1) < SMALL_NUMBER) {
            tan_1 = ZERO_F;
          } else {
            tan_1 = ONE_F/tan_1;
          }

          Pnor_parC = 3.0*etaBragC*(bbgradVC - divVC/3.0);

          tau_xyC = Pnor_parC*biC*bjC;
          tau_zzC = Pnor_parC*(bkC*bkC - 1./3.);

          EXPAND( viscSrc(IDIR,k,j,i) = ZERO_F;  ,
                  viscSrc(JDIR,k,j,i) = (tau_xyC - tau_zzC*tan_1)/x1(i);  ,
                  viscSrc(KDIR,k,j,i) = ZERO_F;  )
        #endif


        // Update flux with the stress tensor
        EXPAND( Flux(MX1, k, j, i) -= tau_xy; ,
                Flux(MX2, k, j, i) -= tau_yy; ,
                Flux(MX3, k, j, i) -= tau_yz; )

        #if HAVE_ENERGY
        Flux(ENG, k, j, i) -= EXPAND(   0.5*(Vc(VX1,k,j,i) + Vc(VX1,k,j-1,i))*tau_xy  ,
                                      + 0.5*(Vc(VX2,k,j,i) + Vc(VX2,k,j-1,i))*tau_yy  ,
                                      + 0.5*(Vc(VX3,k,j,i) + Vc(VX3,k,j-1,i))*tau_yz);
        #endif

        real locdmax = etaBrag/(0.5*(Vc(RHO,k,j,i)+Vc(RHO,k,j-1,i)));
        dMax(k,j,i) = FMAX(dMax(k,j,i),locdmax);
      }


      ///////////////////////////////////////////
      // KDIR sweep                            //
      ///////////////////////////////////////////
      if(dir == KDIR) {
        if(haveViscosity == UserDefFunction) {
          etaBragC = etaBragArr(k,j,i);
          if(haveSlopeLimiter) {
            etaBrag = 2.*(etaBragArr(k-1,j,i)*etaBragArr(k,j,i))/(etaBragArr(k-1,j,i)+etaBragArr(k,j,i));
          } else {
          etaBrag = HALF_F*(etaBragArr(k - 1,j,i)+etaBragArr(k,j,i));
          }
        } else {
          etaBragC = etaBrag = etaBragConstant;
        }

        if (haveSlopeLimiter) {
          dVxi = slopeLimiter(SL_DX(Vc,VX1,k,j,i + 1)/dx1(i+1), SL_DX(Vc,VX1,k - 1,j,i + 1)/dx1(i+1));
          dVxi = slopeLimiter(dVxi, slopeLimiter(SL_DX(Vc,VX1,k,j,i)/dx1(i), SL_DX(Vc,VX1,k - 1,j,i)/dx1(i)));
          dVyi = slopeLimiter(SL_DX(Vc,VX2,k,j,i + 1)/dx1(i+1), SL_DX(Vc,VX2,k - 1,j,i + 1)/dx1(i+1));
          dVyi = slopeLimiter(dVyi, slopeLimiter(SL_DX(Vc,VX2,k,j,i)/dx1(i), SL_DX(Vc,VX2,k - 1,j,i)/dx1(i)));
          dVzi = slopeLimiter(SL_DX(Vc,VX3,k,j,i + 1)/dx1(i+1), SL_DX(Vc,VX3,k - 1,j,i + 1)/dx1(i+1));
          dVzi = slopeLimiter(dVzi, slopeLimiter(SL_DX(Vc,VX3,k,j,i)/dx1(i), SL_DX(Vc,VX3,k - 1,j,i)/dx1(i)));
        } else {
        EXPAND(  dVxi = D_DX_K(Vc,VX1)/dx1(i); ,
                 dVyi = D_DX_K(Vc,VX2)/dx1(i); ,
                 dVzi = D_DX_K(Vc,VX3)/dx1(i); )
        }

        bi = BX_K;
        Bmag = BX_K*BX_K;
        #if DIMENSIONS >= 2
          bj = BY_K;
          Bmag += BY_K*BY_K;
          if (haveSlopeLimiter) {
            dVxj = slopeLimiter(SL_DY(Vc,VX1,k,j + 1,i)/dx2(j+1), SL_DY(Vc,VX1,k - 1,j + 1,i)/dx2(j+1));
            dVxj = slopeLimiter(dVxj, slopeLimiter(SL_DY(Vc,VX1,k,j,i)/dx2(j), SL_DY(Vc,VX1,k - 1,j,i)/dx2(j)));
            dVyj = slopeLimiter(SL_DY(Vc,VX2,k,j + 1,i)/dx2(j+1), SL_DY(Vc,VX2,k - 1,j + 1,i)/dx2(j+1));
            dVyj = slopeLimiter(dVyj, slopeLimiter(SL_DY(Vc,VX2,k,j,i)/dx2(j), SL_DY(Vc,VX2,k - 1,j,i)/dx2(j)));
            dVzj = slopeLimiter(SL_DY(Vc,VX3,k,j + 1,i)/dx2(j+1), SL_DY(Vc,VX3,k - 1,j + 1,i)/dx2(j+1));
            dVzj = slopeLimiter(dVzj, slopeLimiter(SL_DY(Vc,VX3,k,j,i)/dx2(j), SL_DY(Vc,VX3,k - 1,j,i)/dx2(j)));
          } else {
            EXPAND(  dVxj = D_DY_K(Vc,VX1)/dx2(j); ,
                     dVyj = D_DY_K(Vc,VX2)/dx2(j); ,
                     dVzj = D_DY_K(Vc,VX3)/dx2(j); )
          }

          #if DIMENSIONS == 3
              bk = BZ_K;
              Bmag += BZ_K*BZ_K;
              EXPAND (  dVxk = D_DZ_K(Vc,VX1)/dx3(k); ,
                        dVyk = D_DZ_K(Vc,VX2)/dx3(k); ,
                        dVzk = D_DZ_K(Vc,VX3)/dx3(k); )
            }
          #endif
        #endif
        if(Bmag< 0.001*SMALL_NUMBER) {
          Bmag = SQRT(Bmag) + 0.000001*SMALL_NUMBER;
        } else {
          Bmag = SQRT(Bmag);
        }
        bi /= Bmag;
        bj /= Bmag;
        bk /= Bmag;

        #if GEOMETRY == CARTESIAN
          bbgradV = D_EXPAND( bi*bi*dVxi + bj*bi*dVxj + bk*bi*dVxk,
                            + bi*bj*dVyi + bj*bj*dVyj + bk*bj*dVyk,
                            + bi*bk*dVzi + bj*bk*dVzj + bk*bk*dVzk);

          divV = D_EXPAND(dVxi, + dVyj, + dVzk);

          Pnor_par = 3.0*etaBrag*(bbgradV - divV/3.0);

          tau_xz = Pnor_par*bi*bk;
          tau_yz = Pnor_par*bk*bj;
          tau_zz = Pnor_par*(bk*bk - 1./3.);

        // No source term in cartesian geometry
        #endif // GEOMETRY == CARTESIAN


        #if GEOMETRY == SPHERICAL
          real vx1i, vx2i, vx3i;
          vx1i = vx2i = vx3i = ZERO_F;

          vx1i = 0.5*(Vc(VX1,k-1,j,i)+Vc(VX1,k,j,i));
          #if COMPONENTS >= 2
            vx2i = 0.5*(Vc(VX2,k-1,j,i)+Vc(VX2,k,j,i));
            #if COMPONENTS == 3
              vx3i = 0.5*(Vc(VX3,k-1,j,i)+Vc(VX3,k,j,i));
            #endif
          #endif

          bbgradV = D_EXPAND(
      bi*bi*dVxi + bj*bi*1./x1(i)*dVxj + bk*bi*s_1/x1(i)*dVxk,
    + bi*bj*(dVyi - vx2i/x1(i)) + bj*bj*(1./x1(i)*dVyj + vx1i/x1(i)) + bk*bj*s_1/x1(i)*dVyk,
    + bi*bk*(dVzi - vx3i/x1(i)) + bj*bk*(1./x1(i)*dVzj - tan_1/x1(i)*vx3i)
                                        + bk*bk*(s_1/x1(i)*dVzk + vx1i/x1(i) + tan_1/x1(i)*vx2i));

          divV = 2.0*vx1i/x1(i) + dVxi + dVyj/x1(i) + tan_1*vx2i/x1(i) + dVzk/x1(i)*s_1;


          Pnor_par = 3.0*etaBrag*(bbgradV - divV/3.0);

          tau_xz = Pnor_par*bi*bk;
          tau_yz = Pnor_par*bk*bj;
          tau_zz = Pnor_par*(bk*bk - 1./3.);
        #endif

        // Update flux with the stress tensor
        EXPAND( Flux(MX1, k, j, i) -= tau_xz; ,
                Flux(MX2, k, j, i) -= tau_yz; ,
                Flux(MX3, k, j, i) -= tau_zz; )

        #if HAVE_ENERGY
        Flux(ENG, k, j, i) -= EXPAND(   0.5*(Vc(VX1,k,j,i) + Vc(VX1,k-1,j,i))*tau_xz  ,
                                      + 0.5*(Vc(VX2,k,j,i) + Vc(VX2,k-1,j,i))*tau_yz  ,
                                      + 0.5*(Vc(VX3,k,j,i) + Vc(VX3,k-1,j,i))*tau_zz);
        #endif

        real locdmax = etaBrag/(0.5*(Vc(RHO,k,j,i)+Vc(RHO,k-1,j,i)));
        dMax(k,j,i) = FMAX(dMax(k,j,i),locdmax);
      }
  });
  idfx::popRegion();
}
}
