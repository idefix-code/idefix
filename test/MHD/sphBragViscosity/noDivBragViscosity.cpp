// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// This source code is largely inspired from the viscous_flux of Pluto4.2
// ((c) P. Tzeferacos & A. Mignone)

// Implementation of monotonicity-preserving viscous flux following ZuHone et al.,
// ApJ

#include <string>

#include "bragViscosity.hpp"
#include "dataBlock.hpp"
#include "fluid.hpp"


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

////or what should be the exact equivalent, but for some reason isn't, maybe ask Geoffroy why:
//#define   BX_I 0.5*(Vc(BX1,k,j,i) + Vc(BX1,k,j,i - 1))
//#define   BY_J 0.5*(Vc(BX2,k,j,i) + Vc(BX2,k,j - 1,i))
//#define   BZ_K 0.5*(Vc(BX3,k,j,i) + Vc(BX3,k - 1,j,i))
//#define   BY_I 0.5*(Vc(BX2,k,j,i) + Vc(BX2,k,j,i - 1))
//#define   BZ_I 0.5*(Vc(BX3,k,j,i) + Vc(BX3,k,j,i - 1))
//#define   BX_J 0.5*(Vc(BX1,k,j,i) + Vc(BX1,k,j - 1,i))
//#define   BZ_J 0.5*(Vc(BX3,k,j,i) + Vc(BX3,k,j - 1,i))
//#define   BX_K 0.5*(Vc(BX1,k,j,i) + Vc(BX1,k - 1,j,i))
//#define   BY_K 0.5*(Vc(BX2,k,j,i) + Vc(BX2,k - 1,j,i))


void BragViscosity::InitArrays() {
  // Allocate and fill arrays when needed
  #if GEOMETRY != CARTESIAN
    one_dmu = IdefixArray1D<real>("BragViscosity_1dmu", data->np_tot[JDIR]);
    IdefixArray1D<real> dmu = one_dmu;
    IdefixArray1D<real> th = data->x[JDIR];
    idefix_for("ViscousInitGeometry",1,data->np_tot[JDIR],
      KOKKOS_LAMBDA(int j) {
        real scrch =  FABS((1.0-cos(th(j)))*(sin(th(j)) >= 0.0 ? 1.0:-1.0)
                     -(1.0-cos(th(j-1))) * (sin(th(j-1)) > 0.0 ? 1.0:-1.0));
        dmu(j) = 1.0/scrch;
      });
  #endif
  bragViscSrc = IdefixArray4D<real>("BragViscosity_source", COMPONENTS, data->np_tot[KDIR],
                                                                data->np_tot[JDIR],
                                                                data->np_tot[IDIR]);
}
void BragViscosity::ShowConfig() {
  if(status.status==Constant) {
    idfx::cout << "Braginskii Viscosity: ENABLED with constant braginskii viscosity etaBrag="
                    << this->etaBrag << " ."<< std::endl;
  } else if (status.status==UserDefFunction) {
    idfx::cout << "Braginskii Viscosity: ENABLED with user-defined braginskii viscosity function."
                   << std::endl;
    if(!bragViscousDiffusivityFunc) {
      IDEFIX_ERROR("No braginskii viscosity function has been enrolled");
    }
  } else {
    IDEFIX_ERROR("Unknown braginskii viscosity mode");
  }
  if(this->status.isExplicit) {
    idfx::cout << "Braginskii Viscosity: uses an explicit time integration." << std::endl;
  } else if(this->status.isRKL) {
    idfx::cout << "Braginskii Viscosity: uses a Runge-Kutta-Legendre time integration."
                << std::endl;
  } else {
    IDEFIX_ERROR("Unknown time integrator for braginskii viscosity.");
  }
  if(haveSlopeLimiter) {
    idfx::cout << "Braginskii Viscosity: uses a slope limiter." << std::endl;
  }
}

void BragViscosity::EnrollBragViscousDiffusivity(DiffusivityFunc myFunc) {
  if(this->status.status < UserDefFunction) {
    IDEFIX_WARNING("Braginskii viscous diffusivity enrollment requires Hydro/BragViscosity "
                 "to be set to userdef in .ini file");
  }
  this->bragViscousDiffusivityFunc = myFunc;
}

// This function computes the viscous flux and stores it in hydro->fluxRiemann
// (this avoids an extra array)
// Associated source terms, present in non-cartesian geometry are also computed
// and stored in this->viscSrc for later use (in calcRhs).
void BragViscosity::AddBragViscousFlux(int dir, const real t, const IdefixArray4D<real> &Flux) {
  idfx::pushRegion("BragViscosity::AddBragViscousFlux");
  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Vs = this->Vs;
  IdefixArray4D<real> bragViscSrc = this->bragViscSrc;
  IdefixArray3D<real> dMax = this->dMax;
  IdefixArray3D<real> etaBragArr = this->etaBragArr;
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

  #if GEOMETRY == SPHERICAL
  IdefixArray1D<real> sinx2 = this->data->sinx2;
  IdefixArray1D<real> tanx2 = this->data->tanx2;
  IdefixArray1D<real> sinx2m = this->data->sinx2m;
  IdefixArray1D<real> tanx2m = this->data->tanx2m;
  #endif

  HydroModuleStatus haveViscosity = this->status.status;

  // Braginskii Viscosity
  // Compute viscosity if needed
  if(haveViscosity == UserDefFunction && dir == IDIR) {
    if(bragViscousDiffusivityFunc) {
      bragViscousDiffusivityFunc(*this->data, t, etaBragArr);
    } else {
      IDEFIX_ERROR("No user-defined Braginskii viscosity function has been enrolled");
    }
  }

  #if GEOMETRY == CYLINDRICAL || GEOMETRY == POLAR
    IDEFIX_WARNING("Braginskii viscosity in cylindrical and polar coordinates "
                   "has been implemented but never tested so far.");
  #endif

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->data->beg[IDIR];
  iend = this->data->end[IDIR];
  jbeg = this->data->beg[JDIR];
  jend = this->data->end[JDIR];
  kbeg = this->data->beg[KDIR];
  kend = this->data->end[KDIR];

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

      bi = bj = bk = ZERO_F;

      tau_xx = tau_xy = tau_xz = ZERO_F;
      tau_yy = tau_yz = ZERO_F;
      tau_zz = ZERO_F;

      dVxi = dVxj = dVxk = ZERO_F;
      dVyi = dVyj = dVyk = ZERO_F;
      dVzi = dVzj = dVzk = ZERO_F;

      real etaBrag;

      //Cell-centered values to compute source terms
      #if GEOMETRY != CARTESIAN
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

        biC = bjC = bkC = ZERO_F;
        tau_xxC = tau_xyC = tau_xzC = ZERO_F;
        tau_yyC = tau_yzC = ZERO_F;
        tau_zzC = ZERO_F;

        dVxiC = dVxjC = dVxkC = ZERO_F;
        dVyiC = dVyjC = dVykC = ZERO_F;
        dVziC = dVzjC = dVzkC = ZERO_F;

        real etaBragC;

        //Compute values at the center of the cells, no slope limiter are used at this stage
        EXPAND(  dVxiC = HALF_F*(Vc(VX1,k,j,i + 1) - Vc(VX1,k,j,i - 1))/dx1(i); ,
                 dVyiC = HALF_F*(Vc(VX2,k,j,i + 1) - Vc(VX2,k,j,i - 1))/dx1(i); ,
                 dVziC = HALF_F*(Vc(VX3,k,j,i + 1) - Vc(VX3,k,j,i - 1))/dx1(i); )
        EXPAND(  dVxjC = HALF_F*(Vc(VX1,k,j + 1,i) - Vc(VX1,k,j - 1,i))/dx2(j); ,
                 dVyjC = HALF_F*(Vc(VX2,k,j + 1,i) - Vc(VX2,k,j - 1,i))/dx2(j); ,
                 dVzjC = HALF_F*(Vc(VX3,k,j + 1,i) - Vc(VX3,k,j - 1,i))/dx2(j); )
        EXPAND (  dVxkC = HALF_F*(Vc(VX1,k + 1,j,i) - Vc(VX1,k - 1,j,i))/dx3(k); ,
                  dVykC = HALF_F*(Vc(VX2,k + 1,j,i) - Vc(VX2,k - 1,j,i))/dx3(k); ,
                  dVzkC = HALF_F*(Vc(VX3,k + 1,j,i) - Vc(VX3,k - 1,j,i))/dx3(k); )

        EXPAND(  biC = Vc(BX1,k,j,i); ,
                 bjC = Vc(BX2,k,j,i); ,
                 bkC = Vc(BX3,k,j,i); )

        BmagC = EXPAND(  Vc(BX1,k,j,i)*Vc(BX1,k,j,i) ,
                       + Vc(BX2,k,j,i)*Vc(BX2,k,j,i) ,
                       + Vc(BX3,k,j,i)*Vc(BX3,k,j,i) );
        if(BmagC< 0.001*SMALL_NUMBER) {
           BmagC = sqrt(BmagC) + 0.000001*SMALL_NUMBER;
        } else {
           BmagC = sqrt(BmagC);
        }
        biC /= BmagC;
        bjC /= BmagC;
        bkC /= BmagC;
      #endif //GEOMETRY != CARTESIAN

      #if GEOMETRY == CYLINDRICAL
        bbgradVC = D_EXPAND(
                     biC*biC*dVxiC + biC*bjC*dVyiC + biC*bkC*(dVziC  - Vc(VX3,k,j,i)/x1(i)),
                     bjC*biC*dVxjC + bjC*bjC*dVyjC + bjC*bkC*dVzjC + bkC*bkC*Vc(VX1,k,j,i)/x1(i),
                     );

        divVC = D_EXPAND(dVxiC + Vc(VX1,k,j,i)/x1(i) ,
                          + dVyjC ,
                           );
        // No cylindrical geometry in 3D!
      #elif GEOMETRY == POLAR
        bbgradVC = EXPAND(
                    biC*biC*dVxiC                         + bjC*biC*1./x1(i)*dVxjC
                                                                                + bkC*biC*dVxkC,
                  + biC*bjC*(dVyiC - Vc(VX2,k,j,i)/x1(i)) + bjC*bjC*(1./x1(i)*dVyjC + Vc(VX1,k,j,i)/x1(i))
                                                                                + bkC*bjC*dVykC,
                  + biC*bkC*dVziC                         + bjC*bkC*1./x1(i)*dVzjC
                                                                                + bkC*bkC*dVzkC,
                  );

        divVC = D_EXPAND(dVxiC + Vc(VX1,k,j,i)/x1(i) ,
                         + 1./x1(i)*dVyjC ,
                         + dVzkC);
      #elif GEOMETRY == SPHERICAL
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

        bbgradVC = EXPAND(
    biC*biC*dVxiC                         + bjC*biC*1./x1(i)*dVxjC
                                                                          + bkC*biC*s_1/x1(i)*dVxkC,
  + biC*bjC*(dVyiC - Vc(VX2,k,j,i)/x1(i)) + bjC*bjC*(1./x1(i)*dVyjC + Vc(VX1,k,j,i)/x1(i))
                                                                          + bkC*bjC*s_1/x1(i)*dVykC,
  + biC*bkC*(dVziC - Vc(VX3,k,j,i)/x1(i)) + bjC*bkC*(1./x1(i)*dVzjC - tan_1/x1(i)*Vc(VX3,k,j,i))
                    + bkC*bkC*(s_1/x1(i)*dVzkC + Vc(VX1,k,j,i)/x1(i) + tan_1/x1(i)*Vc(VX2,k,j,i)));

        divVC = EXPAND(2.*Vc(VX1,k,j,i)/x1(i) + dVxiC,
                        + dVyjC/x1(i) + tan_1*Vc(VX2,k,j,i)/x1(i),
                        + dVzkC/x1(i)*s_1 );
        divVC = 0.;
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
            #if DIMENSIONS == 3
              dVzj = slopeLimiter(SL_DY(Vc,VX3,k,j + 1,i)/dx2(j+1), SL_DY(Vc,VX3,k,j + 1,i - 1)/dx2(j+1));
              dVzj = slopeLimiter(dVzj, slopeLimiter(SL_DY(Vc,VX3,k,j,i)/dx2(j), SL_DY(Vc,VX3,k,j,i - 1)/dx2(j)));
            #endif
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
          Bmag = sqrt(Bmag) + 0.000001*SMALL_NUMBER;
        } else {
          Bmag = sqrt(Bmag);
        }
        bi /= Bmag;
        bj /= Bmag;
        bk /= Bmag;

        #if GEOMETRY == CARTESIAN
          bbgradV = EXPAND( bi*bi*dVxi + bj*bi*dVxj + bk*bi*dVxk,
                            + bi*bj*dVyi + bj*bj*dVyj + bk*bj*dVyk,
                            + bi*bk*dVzi + bj*bk*dVzj + bk*bk*dVzk);

          divV = EXPAND(dVxi, + dVyj, + dVzk);
        // No source term in cartesian geometry
        #else
          real vx1i, vx2i, vx3i;
          vx1i = vx2i = vx3i = ZERO_F;

          EXPAND(  vx1i = 0.5*(Vc(VX1,k,j,i-1)+Vc(VX1,k,j,i)); ,
                   vx2i = 0.5*(Vc(VX2,k,j,i-1)+Vc(VX2,k,j,i)); ,
                   vx3i = 0.5*(Vc(VX3,k,j,i-1)+Vc(VX3,k,j,i)); )

          Pnor_parC = 3.0*etaBragC*(bbgradVC - divVC/3.0);
        #endif // GEOMETRY != CARTESIAN

        #if GEOMETRY == CYLINDRICAL
          bbgradV = D_EXPAND(
                       bi*bi*dVxi + bi*bj*dVyi + bi*bk*(dVzi  - vx3i/x1l(i)),
                       bj*bi*dVxj + bj*bj*dVyj + bj*bk*dVzj + bk*bk*vx1i/x1l(i),
                       );

          divV = D_EXPAND(dVxi + vx1i/x1l(i) ,
                            + dVyj ,
                             );
          // No cylindrical geometry in 3D!

          //cell-centered values for source terms
          tau_zzC = Pnor_parC*(bkC*bkC - 1./3.);

          EXPAND( bragViscSrc(IDIR,k,j,i) =  -tau_zzC/x1(i);  ,
                  bragViscSrc(JDIR,k,j,i) = ZERO_F;          ,
                  bragViscSrc(KDIR,k,j,i) = ZERO_F;          )

        #elif GEOMETRY == POLAR
          bbgradV = EXPAND(
                      bi*bi*dVxi                         + bj*bi*1./x1l(i)*dVxj
                                                                                  + bk*bi*s_1/x(i)*dVxk,
                    + bi*bj*(dVyi - vx2i/x1l(i)) + bj*bj*(1./x1l(i)*dVyj + vx1i/x1l(i))
                                                                                  + bk*bj*dVyk,
                    + bi*bk*dVzi                         + bj*bk*1./x1l(i)*dVzj
                                                                                  + bk*bk*dVzk,
                    );

          divV = D_EXPAND(dVxi + vx1i/x1l(i) ,
                           + 1./x1l(i)*dVyj ,
                           + dVzk);

          //cell-centered values for source terms
          tau_yyC = Pnor_parC*(bjC*bjC - 1./3.);

          EXPAND( bragViscSrc(IDIR,k,j,i) =  -tau_yyC/x1(i);  ,
                  bragViscSrc(JDIR,k,j,i) = ZERO_F;          ,
                  bragViscSrc(KDIR,k,j,i) = ZERO_F;          )

        #elif GEOMETRY == SPHERICAL
          // NOLINT
          bbgradV = EXPAND(
      bi*bi*dVxi                 + bj*bi*1./x1l(i)*dVxj                 + bk*bi*s_1/x1l(i)*dVxk,
    + bi*bj*(dVyi - vx2i/x1l(i)) + bj*bj*(1./x1l(i)*dVyj + vx1i/x1l(i)) + bk*bj*s_1/x1l(i)*dVyk,
    + bi*bk*(dVzi - vx3i/x1l(i)) + bj*bk*(1./x1l(i)*dVzj - tan_1/x1l(i)*vx3i)
                                      + bk*bk*(s_1/x1l(i)*dVzk + vx1i/x1l(i) + tan_1/x1l(i)*vx2i));

          divV = EXPAND(2.0*vx1i/x1l(i) + dVxi,
                          + dVyj/x1l(i) + tan_1*vx2i/x1l(i),
                          + dVzk/x1l(i)*s_1 );
          divV = 0.;

          //cell-centered values for source terms
          tau_yyC = Pnor_parC*(bjC*bjC - 1./3.);
          tau_zzC = Pnor_parC*(bkC*bkC - 1./3.);

          EXPAND( bragViscSrc(IDIR,k,j,i) =  -(tau_yyC + tau_zzC)/x1(i);  ,
                  bragViscSrc(JDIR,k,j,i) = ZERO_F;          ,
                  bragViscSrc(KDIR,k,j,i) = ZERO_F;          )
        #endif

        Pnor_par = 3.0*etaBrag*(bbgradV - divV/3.0);

        tau_xx = Pnor_par*(bi*bi-1./3.);
        tau_xy = Pnor_par*bi*bj;
        tau_xz = Pnor_par*bi*bk;

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
          #if DIMENSIONS >= 2
            dVyi = slopeLimiter(SL_DX(Vc,VX2,k,j,i + 1)/dx1(i+1), SL_DX(Vc,VX2,k,j - 1,i + 1)/dx1(i+1));
            dVyi = slopeLimiter(dVyi, slopeLimiter(SL_DX(Vc,VX2,k,j,i)/dx1(i), SL_DX(Vc,VX2,k,j - 1,i)/dx1(i)));
            #if DIMENSIONS == 3
              dVzi = slopeLimiter(SL_DX(Vc,VX3,k,j,i + 1)/dx1(i+1), SL_DX(Vc,VX3,k,j - 1,i + 1)/dx1(i+1));
              dVzi = slopeLimiter(dVzi, slopeLimiter(SL_DX(Vc,VX3,k,j,i)/dx1(i), SL_DX(Vc,VX3,k,j - 1,i)/dx1(i)));
            #endif
          #endif
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
          Bmag = sqrt(Bmag) + 0.000001*SMALL_NUMBER;
        } else {
          Bmag = sqrt(Bmag);
        }
        bi /= Bmag;
        bj /= Bmag;
        bk /= Bmag;

        #if GEOMETRY == CARTESIAN
          bbgradV = EXPAND( bi*bi*dVxi + bj*bi*dVxj + bk*bi*dVxk,
                            + bi*bj*dVyi + bj*bj*dVyj + bk*bj*dVyk,
                            + bi*bk*dVzi + bj*bk*dVzj + bk*bk*dVzk);

          divV = EXPAND(dVxi, + dVyj, + dVzk);

          Pnor_par = 3.0*etaBrag*(bbgradV - divV/3.0);

          tau_xy = Pnor_par*bi*bj;
          tau_yy = Pnor_par*(bj*bj - 1./3.);
          tau_yz = Pnor_par*bj*bk;
        // No source term in cartesian geometry
        #else
          real vx1i, vx2i, vx3i;
          vx1i = vx2i = vx3i = ZERO_F;

          EXPAND(  vx1i = 0.5*(Vc(VX1,k,j-1,i)+Vc(VX1,k,j,i)); ,
                   vx2i = 0.5*(Vc(VX2,k,j-1,i)+Vc(VX2,k,j,i)); ,
                   vx3i = 0.5*(Vc(VX3,k,j-1,i)+Vc(VX3,k,j,i)); )

          Pnor_parC = 3.0*etaBragC*(bbgradVC - divVC/3.0);
        #endif // GEOMETRY != CARTESIAN

        #if GEOMETRY == CYLINDRICAL
          bbgradV = D_EXPAND(
                       bi*bi*dVxi + bi*bj*dVyi + bi*bk*(dVzi  - vx3i/x1(i)),
                       bj*bi*dVxj + bj*bj*dVyj + bj*bk*dVzj + bk*bk*vx1i/x1(i),
                       );

          divV = D_EXPAND(dVxi + vx1i/x1(i) ,
                            + dVyj ,
                             );
          // No cylindrical geometry in 3D!
          // No source term along the z-axis in cylindrical/polar coordinates
          EXPAND( bragViscSrc(IDIR,k,j,i) = ZERO_F;  ,
                  bragViscSrc(JDIR,k,j,i) = ZERO_F;  ,
                  bragViscSrc(KDIR,k,j,i) = ZERO_F;  )
        #elif GEOMETRY == POLAR
          bbgradV = EXPAND(
                      bi*bi*dVxi                         + bj*bi*1./x1(i)*dVxj
                                                                                  + bk*bi*s_1/x(i)*dVxk,
                    + bi*bj*(dVyi - vx2i/x1(i)) + bj*bj*(1./x1(i)*dVyj + vx1i/x1(i))
                                                                                  + bk*bj*dVyk,
                    + bi*bk*dVzi                         + bj*bk*1./x1(i)*dVzj
                                                                                  + bk*bk*dVzk,
                    );

          divV = D_EXPAND(dVxi + vx1i/x1(i) ,
                           + 1./x1(i)*dVyj ,
                           + dVzk);

          // See Geoffroy's trick for curvature source terms in cylindrical/polar coordinates
          EXPAND( bragViscSrc(IDIR,k,j,i) = ZERO_F;  ,
                  bragViscSrc(JDIR,k,j,i) = ZERO_F;  ,
                  bragViscSrc(KDIR,k,j,i) = ZERO_F;  )
//          //cell-centered values for source terms
//          tau_xyC = Pnor_parC*biC*bjC;
//
//          EXPAND( bragViscSrc(IDIR,k,j,i) = ZERO_F;  ,
//                  bragViscSrc(JDIR,k,j,i) = tau_xyC/x1(i);  ,
//                  bragViscSrc(KDIR,k,j,i) = ZERO_F;  )

        #elif GEOMETRY == SPHERICAL
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

          bbgradV = EXPAND(
      bi*bi*dVxi + bj*bi*1./x1(i)*dVxj + bk*bi*s_1/x1(i)*dVxk,
    + bi*bj*(dVyi - vx2i/x1(i)) + bj*bj*(1./x1(i)*dVyj + vx1i/x1(i)) + bk*bj*s_1/x1(i)*dVyk,
    + bi*bk*(dVzi - vx3i/x1(i)) + bj*bk*(1./x1(i)*dVzj - tan_1/x1(i)*vx3i)
                                        + bk*bk*(s_1/x1(i)*dVzk + vx1i/x1(i) + tan_1/x1(i)*vx2i));

          divV = EXPAND( 2.0*vx1i/x1(i) + dVxi,
                          +(SIN(x2(j))*Vc(VX2,k,j,i) - FABS(SIN(x2(j-1)))*Vc(VX2,k,j-1,i))/x1(i)
                           *one_dmu(j) ,
                          + dVzk/x1(i)*s_1 );
          divV = 0.;

          //Cell-centered values for the source terms
          tan_1 = tanx2(j);

          // Trick to ensure that the axis does not lead to Nans
          if(FABS(tan_1) < SMALL_NUMBER) {
            tan_1 = ZERO_F;
          } else {
            tan_1 = ONE_F/tan_1;
          }

          tau_xyC = Pnor_parC*biC*bjC;
          tau_zzC = Pnor_parC*(bkC*bkC - 1./3.);

          EXPAND( bragViscSrc(IDIR,k,j,i) = ZERO_F;  ,
                  bragViscSrc(JDIR,k,j,i) = (tau_xyC - tau_zzC*tan_1)/x1(i);  ,
                  bragViscSrc(KDIR,k,j,i) = ZERO_F;  )
        #endif

        Pnor_par = 3.0*etaBrag*(bbgradV - divV/3.0);

        tau_xy = Pnor_par*bi*bj;
        tau_yy = Pnor_par*(bj*bj - 1./3.);
        tau_yz = Pnor_par*bj*bk;

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
          #endif
        #endif
        if(Bmag< 0.001*SMALL_NUMBER) {
          Bmag = sqrt(Bmag) + 0.000001*SMALL_NUMBER;
        } else {
          Bmag = sqrt(Bmag);
        }
        bi /= Bmag;
        bj /= Bmag;
        bk /= Bmag;

        #if GEOMETRY == CARTESIAN
          bbgradV = EXPAND( bi*bi*dVxi + bj*bi*dVxj + bk*bi*dVxk,
                            + bi*bj*dVyi + bj*bj*dVyj + bk*bj*dVyk,
                            + bi*bk*dVzi + bj*bk*dVzj + bk*bk*dVzk);

          divV = EXPAND(dVxi, + dVyj, + dVzk);

          Pnor_par = 3.0*etaBrag*(bbgradV - divV/3.0);

          tau_xz = Pnor_par*bi*bk;
          tau_yz = Pnor_par*bk*bj;
          tau_zz = Pnor_par*(bk*bk - 1./3.);

        // No source term in cartesian geometry
        #else
          real vx1i, vx2i, vx3i;
          vx1i = vx2i = vx3i = ZERO_F;

          EXPAND(  vx1i = 0.5*(Vc(VX1,k-1,j,i)+Vc(VX1,k,j,i)); ,
                   vx2i = 0.5*(Vc(VX2,k-1,j,i)+Vc(VX2,k,j,i)); ,
                   vx3i = 0.5*(Vc(VX3,k-1,j,i)+Vc(VX3,k,j,i)); )
        #endif // GEOMETRY != CARTESIAN

        #if GEOMETRY == CYLINDRICAL
          bbgradV = D_EXPAND(
                       bi*bi*dVxi + bi*bj*dVyi + bi*bk*(dVzi  - vx3i/x1(i)),
                       bj*bi*dVxj + bj*bj*dVyj + bj*bk*dVzj + bk*bk*vx1i/x1(i),
                       );

          divV = D_EXPAND(dVxi + vx1i/x1(i) ,
                            + dVyj ,
                             );
          // No cylindrical geometry in 3D!
          // See Geoffroy's trick for curvature source terms in cylindrical/polar coordinates
          bragViscSrc(IDIR,k,j,i) = ZERO_F;
          bragViscSrc(JDIR,k,j,i) = ZERO_F;
          bragViscSrc(KDIR,k,j,i) = ZERO_F;
        #elif GEOMETRY == POLAR
          bbgradV = EXPAND(
                      bi*bi*dVxi                         + bj*bi*1./x1(i)*dVxj
                                                                                  + bk*bi*s_1/x(i)*dVxk,
                    + bi*bj*(dVyi - vx2i/x1(i)) + bj*bj*(1./x1(i)*dVyj + vx1i/x1(i))
                                                                                  + bk*bj*dVyk,
                    + bi*bk*dVzi                         + bj*bk*1./x1(i)*dVzj
                                                                                  + bk*bk*dVzk,
                    );

          divV = D_EXPAND(dVxi + vx1i/x1(i) ,
                           + 1./x1(i)*dVyj ,
                           + dVzk);

          // No source term along the z-axis in cylindrical/polar coordinates
          bragViscSrc(IDIR,k,j,i) = ZERO_F;
          bragViscSrc(JDIR,k,j,i) = ZERO_F;
          bragViscSrc(KDIR,k,j,i) = ZERO_F;

        #elif GEOMETRY == SPHERICAL
          bbgradV = EXPAND(
      bi*bi*dVxi + bj*bi*1./x1(i)*dVxj + bk*bi*s_1/x1(i)*dVxk,
    + bi*bj*(dVyi - vx2i/x1(i)) + bj*bj*(1./x1(i)*dVyj + vx1i/x1(i)) + bk*bj*s_1/x1(i)*dVyk,
    + bi*bk*(dVzi - vx3i/x1(i)) + bj*bk*(1./x1(i)*dVzj - tan_1/x1(i)*vx3i)
                                        + bk*bk*(s_1/x1(i)*dVzk + vx1i/x1(i) + tan_1/x1(i)*vx2i));

          divV = 2.0*vx1i/x1(i) + dVxi + dVyj/x1(i) + tan_1*vx2i/x1(i) + dVzk/x1(i)*s_1;
          divV = 0.;

          bragViscSrc(IDIR,k,j,i) = ZERO_F;
          bragViscSrc(JDIR,k,j,i) = ZERO_F;
          bragViscSrc(KDIR,k,j,i) = ZERO_F;
        #endif

        Pnor_par = 3.0*etaBrag*(bbgradV - divV/3.0);

        tau_xz = Pnor_par*bi*bk;
        tau_yz = Pnor_par*bk*bj;
        tau_zz = Pnor_par*(bk*bk - 1./3.);

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

real minmodV(const real dvp, const real dvm) {
  real dq= 0.0;
  if(dvp*dvm >0.0) {
    real dq = ( fabs(dvp) < fabs(dvm) ? dvp : dvm);
  }
  return(dq);
}

real monotonizedCentralV(const real dvp, const real dvm) {
  real dq = 0;
  if(dvp*dvm >0.0) {
    real dqc = 0.5*(dvp+dvm);
    real d2q = 2.0*( fabs(dvp) < fabs(dvm) ? dvp : dvm);
    dq= fabs(d2q) < fabs(dqc) ? d2q : dqc;
  }
  return(dq);
}

real vanLeerV(const real dvp, const real dvm) {
  real dq= 0.0;
  dq = (dvp*dvm > ZERO_F ? TWO_F*dvp*dvm/(dvp + dvm) : ZERO_F);
  return(dq);
}
