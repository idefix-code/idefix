// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// This source code is largely inspired from the th_flux of Pluto4.2
// ((c) P. Tzeferacos & A. Mignone)

// Implementation of monotonicity-preserving heat flux following Sharma & Hammett 2007,
// Jour. of Comp. Physics


#include <string>

#include "bragThermalDiffusion.hpp"
#include "dataBlock.hpp"
#include "fluid.hpp"

//ADAPTED TO GET TEMPERATURE DERIVATIVES
#define D_DX_I_tmp(q,n,m)  (q(n,k,j,i)/q(m,k,j,i) - q(n,k,j,i - 1)/q(m,k,j,i - 1))
#define D_DY_J_tmp(q,n,m)  (q(n,k,j,i)/q(m,k,j,i) - q(n,k,j - 1,i)/q(m,k,j - 1,i))
#define D_DZ_K_tmp(q,n,m)  (q(n,k,j,i)/q(m,k,j,i) - q(n,k - 1,j,i)/q(m,k - 1,j,i))

#define SL_DX_tmp(q,n,m,iz,iy,ix)  (q(n,iz,iy,ix)/q(m,iz,iy,ix) - q(n,iz,iy,ix - 1)/q(m,iz,iy,ix - 1))
#define SL_DY_tmp(q,n,m,iz,iy,ix)  (q(n,iz,iy,ix)/q(m,iz,iy,ix) - q(n,iz,iy - 1,ix)/q(m,iz,iy - 1,ix))
#define SL_DZ_tmp(q,n,m,iz,iy,ix)  (q(n,iz,iy,ix)/q(m,iz,iy,ix) - q(n,iz - 1,iy,ix)/q(m,iz - 1,iy,ix))

#define D_DY_I_tmp(q,n,m)  (  0.25*(q(n,k,j + 1,i    ) / q(m,k,j + 1,i)          \
                                  + q(n,k,j + 1,i - 1) / q(m,k,j + 1,i - 1))     \
                            - 0.25*(q(n,k,j - 1,i)     / q(m,k,j - 1,i)          \
                                  + q(n,k,j - 1,i - 1) / q(m,k,j - 1,i - 1)))

#define D_DZ_I_tmp(q,n,m)  (  0.25*(q(n,k + 1,j,i)     / q(m,k + 1,j,i)       \
                                  + q(n,k + 1,j,i - 1) / q(m,k + 1,j,i - 1))  \
                            - 0.25*(q(n,k - 1,j,i)     / q(m,k - 1,j,i)       \
                                  + q(n,k - 1,j,i - 1) / q(m,k - 1,j,i - 1)))

#define D_DX_J_tmp(q,n,m)  (  0.25*(q(n,k,j,i + 1)     / q(m,k,j,i + 1)      \
                                  + q(n,k,j - 1,i + 1) / q(m,k,j - 1,i + 1)) \
                            - 0.25*(q(n,k,j,i - 1)     / q(m,k,j,i - 1)      \
                                  + q(n,k,j - 1,i - 1) / q(m,k,j - 1,i - 1)))

#define D_DZ_J_tmp(q,n,m)  (  0.25*(q(n,k + 1,j,i)     / q(m,k + 1,j,i)       \
                                  + q(n,k + 1,j - 1,i) / q(m,k + 1,j - 1,i))  \
                            - 0.25*(q(n,k - 1,j,i)     / q(m,k - 1,j,i)       \
                                  + q(n,k - 1,j - 1,i) / q(m,k - 1,j - 1,i)))

#define D_DX_K_tmp(q,n,m)  (  0.25*(q(n,k,j,i + 1)     / q(m,k,j,i + 1)      \
                                  + q(n,k - 1,j,i + 1) / q(m,k - 1,j,i + 1)) \
                            - 0.25*(q(n,k,j,i - 1)     / q(m,k,j,i - 1)      \
                                  + q(n,k - 1,j,i - 1) / q(m,k - 1,j,i - 1)))

#define D_DY_K_tmp(q,n,m)  (  0.25*(q(n,k,j + 1,i)     / q(m,k,j + 1,i)       \
                                  + q(n,k - 1,j + 1,i) / q(m,k - 1,j + 1,i))  \
                            - 0.25*(q(n,k,j - 1,i)     / q(m,k,j - 1,i)       \
                                  + q(n,k - 1,j - 1,i) / q(m,k - 1,j - 1,i))) \

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

void BragThermalDiffusion::ShowConfig() {
  if(status.status==Constant) {
    idfx::cout << "Braginskii Thermal Diffusion: ENEBLED with constant diffusivity kpar="
                    << this->kpar << " and knor=" << this->knor << " ."<< std::endl;
  } else if (status.status==UserDefFunction) {
    idfx::cout << "Braginskii Thermal Diffusion: ENABLED with user-defined diffusivity function."
                   << std::endl;
    if(!diffusivityFunc) {
      IDEFIX_ERROR("No braginskii thermal diffusion function has been enrolled");
    }
  } else {
    IDEFIX_ERROR("Unknown braginskii thermal diffusion mode");
  }
  if(status.isExplicit) {
    idfx::cout << "Braginskii Thermal Diffusion: uses an explicit time integration." << std::endl;
  } else if(status.isRKL) {
    idfx::cout << "Braginskii Thermal Diffusion: uses a Runge-Kutta-Legendre time integration."
                << std::endl;
  } else {
    IDEFIX_ERROR("Unknown time integrator for braginskii thermal diffusion.");
  }
  if(haveSlopeLimiter) {
    idfx::cout << "Braginskii Thermal Diffusion: uses a slope limiter." << std::endl;
  }
}

void BragThermalDiffusion::EnrollBragThermalDiffusivity(BragDiffusivityFunc myFunc) {
  if(this->status.status != UserDefFunction) {
    IDEFIX_WARNING("Braginskii thermal diffusivity enrollment requires Hydro/BragThermalDiffusion "
                 "to be set to userdef in .ini file");
  }
  this->diffusivityFunc = myFunc;
}

// This function computes the thermal Flux and stores it in hydro->fluxRiemann
// (this avoids an extra array)
void BragThermalDiffusion::AddBragDiffusiveFlux(int dir, const real t, const IdefixArray4D<real> &Flux) {
#if HAVE_ENERGY == 1
  idfx::pushRegion("BragThermalDiffusion::AddBragDiffusiveFlux");

  real gamma_m1 = this->gamma - ONE_F;
  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Vs = this->Vs;
  IdefixArray3D<real> dMax = this->dMax;

  HydroModuleStatus haveThermalDiffusion = this->status.status;

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

  IdefixArray1D<real> dx = this->data->dx[dir];

  #if GEOMETRY == POLAR
    IdefixArray1D<real> x1 = this->data->x[IDIR];
  #endif
  #if GEOMETRY == SPHERICAL
    IdefixArray1D<real> rt   = this->data->rt;
    IdefixArray1D<real> dmu  = this->data->dmu;
    IdefixArray1D<real> dx2 = this->data->dx[JDIR];
  #endif

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
  real knorConstant = this->knor;
  real kparConstant = this->kpar;
  IdefixArray3D<real> knorArr = this->knorArr;
  IdefixArray3D<real> kparArr = this->kparArr;

  if(haveThermalDiffusion == UserDefFunction && dir == IDIR) {
    if(diffusivityFunc) {
      idfx::pushRegion("UserDef::BragThermalDiffusivityFunction");
      diffusivityFunc(*this->data, t, kparArr, knorArr);
      idfx::popRegion();
    } else {
      IDEFIX_ERROR("No user-defined thermal diffusion function has been enrolled");
    }
  }

  idefix_for("BragDiffusiveFlux",kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      real knor, kpar;
      real bgradT, Bmag, bn;
      real dT_mag;
      real q;
      real Bi, Bj, Bk, Bn = 0;

      real dTi, dTj, dTk, dTn = 0;
      real locdmax = 0;

      ///////////////////////////////////////////
      // IDIR sweep                            //
      ///////////////////////////////////////////

      if(dir == IDIR) {
        if(haveThermalDiffusion == UserDefFunction) {
          knor = HALF_F*(knorArr(k,j,i-1)+knorArr(k,j,i));
          if(haveSlopeLimiter) {
            kpar = 2.*(kparArr(k,j,i-1)*kparArr(k,j,i))/(kparArr(k,j,i-1)+kparArr(k,j,i));
          } else {
            kpar = HALF_F*(kparArr(k,j,i-1)+kparArr(k,j,i));
          }
        } else {
          knor = knorConstant;
          kpar = kparConstant;
        }

        #if DIMENSIONS == 3
        EXPAND( Bi = BX_I; ,
                Bj = BY_I; ,
                Bk = BZ_I; )
        #elif DIMENSIONS == 2
        EXPAND( Bi = BX_I; ,
                Bj = BY_I; ,
                Bk = ZERO_F;)
        #elif DIMENSIONS == 1
        EXPAND( Bi = BX_I; ,
                Bj = ZERO_F; ,
                Bk = ZERO_F;)
        #endif

        Bn = BX_I;

        #if GEOMETRY == CARTESIAN
          dTi = D_DX_I_tmp(Vc,PRS,RHO)/dx1(i);
          #if DIMENSIONS >= 2
            if (haveSlopeLimiter) {
              dTj = slopeLimiter(SL_DY_tmp(Vc,PRS,RHO,k,j,i-1)/dx2(j),SL_DY_tmp(Vc,PRS,RHO,k,j+1,i-1)/dx2(j+1));
              dTj = slopeLimiter(dTj, slopeLimiter(SL_DY_tmp(Vc,PRS,RHO,k,j,i)/dx2(j),SL_DY_tmp(Vc,PRS,RHO,k,j+1,i)/dx2(j+1)));
            } else {
              dTj = D_DY_I_tmp(Vc,PRS,RHO)/dx2(j);
            }
            #if DIMENSIONS == 3
              if (haveSlopeLimiter) {
                dTk = slopeLimiter(SL_DZ_tmp(Vc,PRS,RHO,k,j,i-1)/dx3(k),SL_DZ_tmp(Vc,PRS,RHO,k+1,j,i-1)/dx3(k+1));
                dTk = slopeLimiter(dTk, slopeLimiter(SL_DZ_tmp(Vc,PRS,RHO,k,j,i)/dx3(k),SL_DZ_tmp(Vc,PRS,RHO,k+1,j,i)/dx3(k+1)));
              } else {
                dTk = D_DZ_I_tmp(Vc,PRS,RHO)/dx3(k);
              }
            #endif
          #endif
        #elif GEOMETRY == CYLINDRICAL
          dTi = D_DX_I_tmp(Vc,PRS,RHO)/dx1(i);
          #if DIMENSIONS >= 2
            if (haveSlopeLimiter) {
              dTj = slopeLimiter(SL_DY_tmp(Vc,PRS,RHO,k,j,i-1)/dx2(j),SL_DY_tmp(Vc,PRS,RHO,k,j+1,i-1)/dx2(j+1));
              dTj = slopeLimiter(dTj, slopeLimiter(SL_DY_tmp(Vc,PRS,RHO,k,j,i)/dx2(j),SL_DY_tmp(Vc,PRS,RHO,k,j+1,i)/dx2(j+1)));
            } else {
              dTj = D_DY_I_tmp(Vc,PRS,RHO)/dx2(j);
            }
          #endif
          // No cylindrical geometry in 3D!
        #elif GEOMETRY == POLAR
          dTi = D_DX_I_tmp(Vc,PRS,RHO)/dx1(i);
          #if DIMENSIONS >= 2
            if (haveSlopeLimiter) {
              dTj = 1./x1(i-1)*slopeLimiter(SL_DY_tmp(Vc,PRS,RHO,k,j,i-1)/dx2(j),SL_DY_tmp(Vc,PRS,RHO,k,j+1,i-1)/dx2(j+1));
              dTj = slopeLimiter(dTj, 1./x1(i)*slopeLimiter(SL_DY_tmp(Vc,PRS,RHO,k,j,i)/dx2(j),SL_DY_tmp(Vc,PRS,RHO,k,j+1,i)/dx2(j+1)));
            } else {
              dTj = 1./x1l(i)*D_DY_I_tmp(Vc,PRS,RHO)/dx2(j); // 1/r dTj/dxj
            }
            #if DIMENSIONS == 3
              if (haveSlopeLimiter) {
                dTk = slopeLimiter(SL_DZ_tmp(Vc,PRS,RHO,k,j,i-1)/dx3(k),SL_DZ_tmp(Vc,PRS,RHO,k+1,j,i-1)/dx3(k+1));
                dTk = slopeLimiter(dTk, slopeLimiter(SL_DZ_tmp(Vc,PRS,RHO,k,j,i)/dx3(k),SL_DZ_tmp(Vc,PRS,RHO,k+1,j,i)/dx3(k+1)));
              } else {
                dTk = D_DZ_I_tmp(Vc,PRS,RHO)/dx3(k);
              }
            #endif
          #endif
        #elif GEOMETRY == SPHERICAL
            real s_1 = sinx2(j);

            // Trick to ensure that the axis does not lead to Nans
          if(FABS(s_1) < SMALL_NUMBER) {
              s_1 = ZERO_F;
            } else {
              s_1 = ONE_F/s_1;
            }

          dTi = D_DX_I_tmp(Vc,PRS,RHO)/dx1(i);
          #if DIMENSIONS >= 2
            if (haveSlopeLimiter) {
              dTj = 1./x1(i-1)*slopeLimiter(SL_DY_tmp(Vc,PRS,RHO,k,j,i-1)/dx2(j),SL_DY_tmp(Vc,PRS,RHO,k,j+1,i-1)/dx2(j+1));
              dTj = slopeLimiter(dTj, 1./x1(i)*slopeLimiter(SL_DY_tmp(Vc,PRS,RHO,k,j,i)/dx2(j),SL_DY_tmp(Vc,PRS,RHO,k,j+1,i)/dx2(j+1)));
            } else {
              dTj = 1./x1l(i)*D_DY_I_tmp(Vc,PRS,RHO)/dx2(j); // 1/r dTj/dxj
            }
            #if DIMENSIONS == 3
              if (haveSlopeLimiter) {
                dTk = s_1/x1(i-1)*slopeLimiter(SL_DZ_tmp(Vc,PRS,RHO,k,j,i-1)/dx3(k),SL_DZ_tmp(Vc,PRS,RHO,k+1,j,i-1)/dx3(k+1));
                dTk = slopeLimiter(dTk, s_1/x1(i)*slopeLimiter(SL_DZ_tmp(Vc,PRS,RHO,k,j,i)/dx3(k),SL_DZ_tmp(Vc,PRS,RHO,k+1,j,i)/dx3(k+1)));
              } else {
                dTk = s_1/x1l(i)*D_DZ_I_tmp(Vc,PRS,RHO)/dx3(k);
              }
            #endif
          #endif
        #endif // GEOMETRY

        locdmax = FMAX(kpar,knor)/(0.5*(Vc(RHO,k,j,i) + Vc(RHO,k,j,i - 1)))*gamma_m1;

        dTn = dTi;
      } else if(dir == JDIR) {
        //////////////
        // JDIR
        /////////////

        if(haveThermalDiffusion == UserDefFunction) {
          knor = HALF_F*(knorArr(k,j-1,i)+knorArr(k,j,i));
          if(haveSlopeLimiter) {
            kpar = 2.*(kparArr(k,j-1,i)*kparArr(k,j,i))/(kparArr(k,j-1,i)+kparArr(k,j,i));
          } else {
            kpar = HALF_F*(kparArr(k,j-1,i)+kparArr(k,j,i));
          }
        } else {
          knor = knorConstant;
          kpar = kparConstant;
        }

        #if DIMENSIONS == 3
        EXPAND( Bi = BX_J; ,
                Bj = BY_J; ,
                Bk = BZ_J; )
        #elif DIMENSIONS == 2
        EXPAND( Bi = BX_J; ,
                Bj = BY_J; ,
                Bk = ZERO_F;)
        #endif

        Bn = BY_J;

        #if GEOMETRY == CARTESIAN
          if (haveSlopeLimiter) {
            dTi = slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k,j-1,i)/dx1(i),SL_DX_tmp(Vc,PRS,RHO,k,j-1,i+1)/dx1(i+1));
            dTi = slopeLimiter(dTi, slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k,j,i)/dx1(i),SL_DX_tmp(Vc,PRS,RHO,k,j,i+1)/dx1(i+1)));
          } else {
            dTi = D_DX_J_tmp(Vc,PRS,RHO)/dx1(i);
          }
          dTj = D_DY_J_tmp(Vc,PRS,RHO)/dx2(j);
          #if DIMENSIONS == 3
            if (haveSlopeLimiter) {
              dTk = slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k,j-1,i)/dx3(k),SL_DX_tmp(Vc,PRS,RHO,k+1,j-1,i)/dx3(k+1));
              dTk = slopeLimiter(dTk, slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k,j,i)/dx3(k),SL_DX_tmp(Vc,PRS,RHO,k+1,j,i)/dx3(k+1)));
            } else {
              dTk = D_DZ_J_tmp(Vc,PRS,RHO)/dx3(k);
            }
          #endif
        #elif GEOMETRY == CYLINDRICAL
          if (haveSlopeLimiter) {
            dTi = slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k,j-1,i)/dx1(i),SL_DX_tmp(Vc,PRS,RHO,k,j-1,i+1)/dx1(i+1));
            dTi = slopeLimiter(dTi, slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k,j,i)/dx1(i),SL_DX_tmp(Vc,PRS,RHO,k,j,i+1)/dx1(i+1)));
          } else {
            dTi = D_DX_J_tmp(Vc,PRS,RHO)/dx1(i);
          }
          dTj = D_DY_J_tmp(Vc,PRS,RHO)/dx2(j);
          // No cylindrical geometry in 3D!
        #elif GEOMETRY == POLAR
          //gradT = ... + 1/r dT/dtheta + ...
          if (haveSlopeLimiter) {
            dTi = slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k,j-1,i)/dx1(i),SL_DX_tmp(Vc,PRS,RHO,k,j-1,i+1)/dx1(i+1));
            dTi = slopeLimiter(dTi, slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k,j,i)/dx1(i),SL_DX_tmp(Vc,PRS,RHO,k,j,i+1)/dx1(i+1)));
          } else {
            dTi = D_DX_J_tmp(Vc,PRS,RHO)/dx1(i);
          }
          dTj = 1./x1(i)*D_DY_J_tmp(Vc,PRS,RHO)/dx2(j);
          #if DIMENSIONS == 3
            if (haveSlopeLimiter) {
              dTk = slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k,j-1,i)/dx3(k),SL_DX_tmp(Vc,PRS,RHO,k+1,j-1,i)/dx3(k+1));
              dTk = slopeLimiter(dTk, slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k,j,i)/dx3(k),SL_DX_tmp(Vc,PRS,RHO,k+1,j,i)/dx3(k+1)));
            } else {
              dTk = D_DZ_J_tmp(Vc,PRS,RHO)/dx3(k);
#                dTk = D_DZ_J_tmp(0,RHO)/dx3(k);
            }
          #endif
        #elif GEOMETRY == SPHERICAL
          if (haveSlopeLimiter) {
            dTi = slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k,j-1,i)/dx1(i),SL_DX_tmp(Vc,PRS,RHO,k,j-1,i+1)/dx1(i+1));
            dTi = slopeLimiter(dTi, slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k,j,i)/dx1(i),SL_DX_tmp(Vc,PRS,RHO,k,j,i+1)/dx1(i+1)));
          } else {
            dTi = D_DX_J_tmp(Vc,PRS,RHO)/dx1(i);
          }
          dTj = 1./x1(i)*D_DY_J_tmp(Vc,PRS,RHO)/dx2(j);
          #if DIMENSIONS == 3
            if (haveSlopeLimiter) {
              real s_1 = sinx2(j-1);

              // Trick to ensure that the axis does not lead to Nans
              if(FABS(s_1) < SMALL_NUMBER) {
                s_1 = ZERO_F;
              } else {
                s_1 = ONE_F/s_1;
              }
              dTk = s_1/x1(i)*slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k,j-1,i)/dx3(k),SL_DX_tmp(Vc,PRS,RHO,k+1,j-1,i)/dx3(k+1));

              s_1 = sinx2(j);

              // Trick to ensure that the axis does not lead to Nans
              if(FABS(s_1) < SMALL_NUMBER) {
                s_1 = ZERO_F;
              } else {
                s_1 = ONE_F/s_1;
              }
              dTk = s_1/x1(i)*slopeLimiter(dTk, slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k,j,i)/dx3(k),SL_DX_tmp(Vc,PRS,RHO,k+1,j,i)/dx3(k+1)));
            } else {
              real s_1 = sinx2m(j);

              // Trick to ensure that the axis does not lead to Nans
              if(FABS(s_1) < SMALL_NUMBER) {
                s_1 = ZERO_F;
              } else {
                s_1 = ONE_F/s_1;
              }

              dTk = s_1/x1(i)*D_DZ_J_tmp(Vc,PRS,RHO)/dx3(k);
            }
          #endif
        #endif // GEOMETRY

        locdmax = FMAX(kpar,knor)/(0.5*(Vc(RHO,k,j,i) + Vc(RHO,k,j - 1,i)))*gamma_m1;

        dTn = dTj;
      } else if(dir == KDIR) {
      //////////////
      // KDIR
      /////////////
        if(haveThermalDiffusion == UserDefFunction) {
          knor = HALF_F*(knorArr(k-1,j,i)+knorArr(k,j,i));
          if(haveSlopeLimiter) {
            kpar = 2.*(kparArr(k-1,j,i)*kparArr(k,j,i))/(kparArr(k-1,j,i)+kparArr(k,j,i));
          } else {
            kpar = HALF_F*(kparArr(k-1,j,i)+kparArr(k,j,i));
          }
        } else {
          knor = knorConstant;
          kpar = kparConstant;
        }

        Bi = BX_K;
        Bj = BY_K;
        Bk = BZ_K;

        Bn = Bk;

        #if GEOMETRY == CARTESIAN
          if (haveSlopeLimiter) {
            dTi = slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k-1,j,i)/dx1(i),SL_DX_tmp(Vc,PRS,RHO,k-1,j,i+1)/dx1(i+1));
            dTi = slopeLimiter(dTi, slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k,j,i)/dx1(i),SL_DX_tmp(Vc,PRS,RHO,k,j,i+1)/dx1(i+1)));
            dTj = slopeLimiter(SL_DY_tmp(Vc,PRS,RHO,k-1,j,i)/dx2(j),SL_DY_tmp(Vc,PRS,RHO,k-1,j+1,i)/dx2(j+1));
            dTj = slopeLimiter(dTj, slopeLimiter(SL_DY_tmp(Vc,PRS,RHO,k,j,i)/dx2(j),SL_DY_tmp(Vc,PRS,RHO,k,j+1,i)/dx2(j+1)));
          } else {
            dTi = D_DX_K_tmp(Vc,PRS,RHO)/dx1(i);
            dTj = D_DY_K_tmp(Vc,PRS,RHO)/dx2(j);
          }
          dTk = D_DZ_K_tmp(Vc,PRS,RHO)/dx3(k);
        // No cylindrical geometry in 3D!
        #elif GEOMETRY == POLAR
          if (haveSlopeLimiter) {
            dTi = slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k-1,j,i)/dx1(i),SL_DX_tmp(Vc,PRS,RHO,k-1,j,i+1)/dx1(i+1));
            dTi = slopeLimiter(dTi, slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k,j,i)/dx1(i),SL_DX_tmp(Vc,PRS,RHO,k,j,i+1)/dx1(i+1)));
            dTj = slopeLimiter(SL_DY_tmp(Vc,PRS,RHO,k-1,j,i)/dx2(j),SL_DY_tmp(Vc,PRS,RHO,k-1,j+1,i)/dx2(j+1));
            dTj = slopeLimiter(dTj, slopeLimiter(SL_DY_tmp(Vc,PRS,RHO,k,j,i)/dx2(j),SL_DY_tmp(Vc,PRS,RHO,k,j+1,i)/dx2(j+1)));
            dTj *= 1./x1(i);
          } else {
            dTi = D_DX_K_tmp(Vc,PRS,RHO)/dx1(i);
            dTj = 1./x1(i)*D_DY_K_tmp(Vc,PRS,RHO)/dx2(j); // 1/r dTj/dxj
          }
          dTk = D_DZ_K_tmp(Vc,PRS,RHO)/dx3(k);
        #elif GEOMETRY == SPHERICAL
            real s_1 = sinx2(j);

          // Trick to ensure that the axis does not lead to Nans
          if(FABS(s_1) < SMALL_NUMBER) {
            s_1 = ZERO_F;
          } else {
            s_1 = ONE_F/s_1;
          }

          if (haveSlopeLimiter) {
            dTi = slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k-1,j,i)/dx1(i),SL_DX_tmp(Vc,PRS,RHO,k-1,j,i+1)/dx1(i+1));
            dTi = slopeLimiter(dTi, slopeLimiter(SL_DX_tmp(Vc,PRS,RHO,k,j,i)/dx1(i),SL_DX_tmp(Vc,PRS,RHO,k,j,i+1)/dx1(i+1)));
            dTj = slopeLimiter(SL_DY_tmp(Vc,PRS,RHO,k-1,j,i)/dx2(j),SL_DY_tmp(Vc,PRS,RHO,k-1,j+1,i)/dx2(j+1));
            dTj = slopeLimiter(dTj, slopeLimiter(SL_DY_tmp(Vc,PRS,RHO,k,j,i)/dx2(j),SL_DY_tmp(Vc,PRS,RHO,k,j+1,i)/dx2(j+1)));
            dTj *= 1./x1(i);
          } else {
            dTi = D_DX_K_tmp(Vc,PRS,RHO)/dx1(i);
            dTj = 1./x1(i)*D_DY_K_tmp(Vc,PRS,RHO)/dx2(j); // 1/r dTj/dxj
          }
          dTk = s_1/x1(i)*D_DZ_K_tmp(Vc,PRS,RHO)/dx3(k);
          //gradT = ... + ... + 1/(r*sin(theta)) dTphi/dphi
        #endif // GEOMETRY

        locdmax = FMAX(kpar,knor)/(0.5*(Vc(RHO,k,j,i) + Vc(RHO,k-1,j,i)))*gamma_m1;

        dTn = dTk;
      }
      // From here, gradients and normal have been computed, so we just need to get the fluxes

      dT_mag  = D_EXPAND(   dTi*dTi,
                        +   dTj*dTj,
                        +   dTk*dTk);

      dT_mag = sqrt(dT_mag) + 0.01*SMALL_NUMBER;

      #if DIMENSIONS == 3
      bgradT = EXPAND( Bi*dTi , + Bj*dTj, +Bk*dTk);
      Bmag = EXPAND( Bi*Bi , + Bj*Bj, + Bk*Bk);
      #elif DIMENSIONS == 2
      bgradT = EXPAND( Bi*dTi , + Bj*dTj, +0);
      Bmag = EXPAND( Bi*Bi , + Bj*Bj, + 0);
      #elif DIMENSIONS == 1
      bgradT = EXPAND( Bi*dTi , + 0, +0);
      Bmag = EXPAND( Bi*Bi , + 0, + 0);
      #endif
      Bmag = sqrt(Bmag);
      Bmag = FMAX(1e-6*SMALL_NUMBER,Bmag);

      bgradT /= Bmag;

      bn = Bn/Bmag; /* -- unit vector component -- */
      q = kpar*bgradT*bn + knor*(dTn - bn*bgradT);
      //qi_mag = sqrt(  (kpar*kpar - knor*knor)*bgradT*bgradT
      //  + knor*knor*dT_mag*dT_mag);

      Flux(ENG, k, j, i) -= q;

//      dMax(k,j,i) += FMAX(dMax(k,j,i),locdmax);
      dMax(k,j,i) = FMAX(dMax(k,j,i),locdmax);
    });
  idfx::popRegion();
#endif // HAVE_ENERGY
}

real minmodTh(const real dvp, const real dvm) {
  real dq= 0.0;
  if(dvp*dvm >0.0) {
    real dq = ( fabs(dvp) < fabs(dvm) ? dvp : dvm);
  }   
  return(dq);
}

real monotonizedCentralTh(const real dvp, const real dvm) {
  real dq = 0;
  if(dvp*dvm >0.0) {
    real dqc = 0.5*(dvp+dvm);
    real d2q = 2.0*( fabs(dvp) < fabs(dvm) ? dvp : dvm);
    dq= fabs(d2q) < fabs(dqc) ? d2q : dqc;
  }
  return(dq);
}

real vanLeerTh(const real dvp, const real dvm) {
  real dq= 0.0;
  dq = (dvp*dvm > ZERO_F ? TWO_F*dvp*dvm/(dvp + dvm) : ZERO_F);
  return(dq);
}

