// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// This source code is largely inspired from the th_flux of Pluto4.2
// ((c) P. Tzeferacos & A. Mignone)
//
//
// ONGOING SPHERICAL K DIRECTION
// 
// Implementation of monotonicity-preserving heat flux following Sharma & Hammett 2007,
// Jour. of Comp. Physics

//#if HAVE_ENERGY

#include <string>

#include "thConductivity.hpp"
#include "dataBlock.hpp"

//ADAPT TO GET TEMPERATURE DERIVATIVES
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


ThConductivity::ThConductivity() {
  // Default constructor
}

#if HAVE_ENERGY
void ThConductivity::Init(Input &input, Grid &grid, Hydro *hydroin) {
  idfx::pushRegion("thConductivity::Init");
  // Save the parent hydro object
  this->hydro = hydroin;
  this->haveBraginskiiConductivity = false;
  this->haveSlopeLimiter = false;

  if(input.CheckEntry("Hydro","thConductivity")>=1) {
    if(input.Get<std::string>("Hydro","thConductivity",1).compare("constant") == 0) {
      this->kappa = input.Get<real>("Hydro","thConductivity",2);
      this->haveThConductivity = Constant;
    } else if(input.Get<std::string>("Hydro","thConductivity",1).compare("userdef") == 0) {
      this->haveThConductivity = UserDefFunction;
      this->kappaArr = IdefixArray3D<real>("ThconductivityKnorArray",hydro->data->np_tot[KDIR],
                                                               hydro->data->np_tot[JDIR],
                                                               hydro->data->np_tot[IDIR]);
    } else if(input.Get<std::string>("Hydro", "thConductivity",1).compare("braginskii") == 0) {
      #if MHD == NO
        IDEFIX_ERROR("MHD is not enabled."
                     "Can only enabled Braginskii conductivity when MHD is enabled");
      #endif
      this->haveBraginskiiConductivity = true;
      if(input.Get<std::string>("Hydro","thConductivity",2).compare("constant") == 0) {
        this->haveThConductivity = Constant;
        this->kpar = input.Get<real>("Hydro","thConductivity",3);
        this->knor = input.GetOrSet<real>("Hydro","thConductivity",4, 0.0);
      } else if(input.Get<std::string>("Hydro","thConductivity", 2).compare("userdef") == 0) {
        this->haveThConductivity = UserDefFunction;
        this->kparArr = IdefixArray3D<real>("ThconductivityKparArray",hydro->data->np_tot[KDIR],
                                                                hydro->data->np_tot[JDIR],
                                                                hydro->data->np_tot[IDIR]);
        this->knorArr = IdefixArray3D<real>("ThconductivityKnorArray",hydro->data->np_tot[KDIR],
                                                                hydro->data->np_tot[JDIR],
                                                                hydro->data->np_tot[IDIR]);
      } else {
        IDEFIX_ERROR("Unknown Braginskii thermal conductivity definition in idefix.ini. "
                     "Can only be constant or userdef.");
      }
    } else if(input.Get<std::string>("Hydro", "thConductivity",1).compare("braginskiilim") == 0) {
      #if MHD == NO
        IDEFIX_ERROR("MHD is not enabled."
                     "Can only enabled BraginskiiLim conductivity when MHD is enabled");
      #endif
      this->haveBraginskiiConductivity = true;
      this->haveSlopeLimiter = true;
      if(input.Get<std::string>("Hydro", "thConductivity",2).compare("minmod") == 0) {
        this->slopeLimiter = minmodTh;
      }
      else if(input.Get<std::string>("Hydro", "thConductivity",2).compare("vanleer") == 0) {
        this->slopeLimiter = vanLeerTh;
      }
      else if(input.Get<std::string>("Hydro", "thConductivity",2).compare("mclim") == 0) {
        this->slopeLimiter = monotonizedCentralTh;
      } else {
        IDEFIX_ERROR("Unknown BraginskiiLim slope limiter definition in idefix.ini. "
                     "Can only be minmod, vanleer or mclim.");
      }
      if(input.Get<std::string>("Hydro","thConductivity",3).compare("constant") == 0) {
        this->haveThConductivity = Constant;
        this->kpar = input.Get<real>("Hydro","thConductivity",4);
        this->knor = input.GetOrSet<real>("Hydro","thConductivity",5, 0.0);
      } else if(input.Get<std::string>("Hydro","thConductivity", 3).compare("userdef") == 0) {
        this->haveThConductivity = UserDefFunction;
        this->kparArr = IdefixArray3D<real>("ThconductivityKnorArray",hydro->data->np_tot[KDIR],
                                                                hydro->data->np_tot[JDIR],
                                                                hydro->data->np_tot[IDIR]);
        this->knorArr = IdefixArray3D<real>("ThconductivityKnorArray",hydro->data->np_tot[KDIR],
                                                                hydro->data->np_tot[JDIR],
                                                                hydro->data->np_tot[IDIR]);
      } else {
        IDEFIX_ERROR("Unknown BraginskiiLim thermal conductivity definition in idefix.ini. "
                     "Can only be constant or userdef.");
      }
    } else {
      IDEFIX_ERROR("Unknown thermal conductivity definition in idefix.ini. "
                   "Can only be constant, userdef, braginskii or braginskiilim.");
    }
  } else {
    IDEFIX_ERROR("I cannot create a ThConductivity object without thConductivity defined"
                   "in the .ini file");
  }

  idfx::popRegion();
}


void ThConductivity::ShowConfig() {
  if(haveThConductivity==Constant) {
    if(haveBraginskiiConductivity) {
      idfx::cout << "ThConductivity: ENABLED with constant "
                 << "Braginskii conductivity kappaPar=" << this->kpar
                << " and kappaNor=" << this->knor << " ." << std::endl;
    } else {
      idfx::cout << "ThConductivity: ENABLED with constant Braginskii conductivity kappa="
                      << this->kappa << " ."<< std::endl;
    }
  } else if (haveThConductivity==UserDefFunction) {
    if(haveBraginskiiConductivity) {
      idfx::cout << "ThConductivity: ENABLED with user-defined "
                 << "Braginskii conductivity function." << std::endl;
      if(!bragThermalConductivityFunc) {
        IDEFIX_ERROR("No Braginskii thermal diffusion function has been enrolled");
      }
    } else {
      idfx::cout << "ThConductivity: ENABLED with user-defined diffusivity function."
                   << std::endl;
      if(!thermalConductivityFunc) {
        IDEFIX_ERROR("No thermal diffusion function has been enrolled");
      }
    }
  } else {
    IDEFIX_ERROR("Unknown thermal conductivity mode");
  }
  if(hydro->thConductivityStatus.isExplicit) {
    idfx::cout << "ThConductivity: uses an explicit time integration." << std::endl;
  } else if(hydro->thConductivityStatus.isRKL) {
    idfx::cout << "ThConductivity: uses a Runge-Kutta-Legendre time integration."
                << std::endl;
  } else {
    IDEFIX_ERROR("Unknown time integrator for thermal conductivity.");
  }
  if(haveSlopeLimiter) {
    idfx::cout << "ThConductivity: uses a slope limiter." << std::endl;
  }
}

void ThConductivity::EnrollThermalConductivity(ThermalConductivityFunc myFunc) {
  if(this->haveThConductivity < UserDefFunction) {
    IDEFIX_ERROR("Thermal conductivity enrollment requires Hydro/ThConductivity "
                 "to be set to userdef in .ini file");
  }
  if(haveBraginskiiConductivity) {
    IDEFIX_ERROR("Use EnrollBragThermalConductivity for Braginskii thermal conduction");
  }
  this->thermalConductivityFunc = myFunc;
}

void ThConductivity::EnrollBragThermalConductivity(BragThermalConductivityFunc myFunc) {
  if(this->haveThConductivity < UserDefFunction) {
    IDEFIX_ERROR("Thermal conductivity enrollment requires Hydro/ThConductivity "
                 "to be set to userdef in .ini file");
  }
  if(!haveBraginskiiConductivity) {
    IDEFIX_ERROR("Use EnrollThermalConductivity for isotropic thermal conduction");
  }
  this->bragThermalConductivityFunc = myFunc;
}

// This function computes the thermal Flux and stores it in hydro->fluxRiemann
// (this avoids an extra array)
void ThConductivity::AddThermalFlux(int dir, const real t) {
  idfx::pushRegion("ThConductivity::AddThermalFlux");

  real gamma_m1 = this->hydro->gamma - ONE_F;
  IdefixArray4D<real> Vc = this->hydro->Vc;
  IdefixArray4D<real> Vs = this->hydro->Vs;
  IdefixArray4D<real> Flux = this->hydro->FluxRiemann;
  IdefixArray3D<real> dMax = this->hydro->dMax;

  HydroModuleStatus haveThConductivity = this->haveThConductivity;

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

//Isotropic case
if (! this->haveBraginskiiConductivity) {
  IdefixArray3D<real> kappaArr = this->kappaArr;
  IdefixArray1D<real> dx = this->hydro->data->dx[dir];

  #if GEOMETRY == POLAR
    IdefixArray1D<real> x1 = this->hydro->data->x[IDIR];
  #endif
  #if GEOMETRY == SPHERICAL
    IdefixArray1D<real> rt   = this->hydro->data->rt;
    IdefixArray1D<real> dmu  = this->hydro->data->dmu;
    IdefixArray1D<real> dx2 = this->hydro->data->dx[JDIR];
  #endif
    // Compute thermal diffusion if needed
  if(haveThConductivity == UserDefFunction && dir == IDIR) {
    if(thermalConductivityFunc) {
      idfx::pushRegion("UserDef::ThermalDiffusivityFunction");
      thermalConductivityFunc(*this->hydro->data, t, kappaArr);
      idfx::popRegion();
    } else {
      IDEFIX_ERROR("No user-defined thermal diffusion function has been enrolled");
    }
  }
  real kappaConstant = this->kappa;
  const int ioffset = (dir==IDIR) ? 1 : 0;
  const int joffset = (dir==JDIR) ? 1 : 0;
  const int koffset = (dir==KDIR) ? 1 : 0;

  idefix_for("ThermalDiffusionFlux",kbeg, kend, jbeg, jend, ibeg, iend,
      KOKKOS_LAMBDA (int k, int j, int i) {
        // Compute gradT
        real gradT;

        gradT = Vc(PRS,k,j,i) / Vc(RHO,k,j,i)
               - Vc(PRS,k-koffset,j-joffset,i-ioffset) / Vc(RHO,k-koffset,j-joffset,i-ioffset);

        // index along dir
        const int ig = ioffset*i + joffset*j + koffset*k;

        // dx at the interface is the averaged between the two adjacent centered dx
        real dl = HALF_F*(dx(ig-1) + dx(ig));
        #if GEOMETRY == POLAR
        if(dir==JDIR)
          dl = dl*x1(i);

        #elif GEOMETRY == SPHERICAL
          if(dir==JDIR)
            dl = dl*rt(i);
          else
            if(dir==KDIR)
              dl = dl*rt(i)*dmu(j)/dx2(j);
        #endif // GEOMETRY

        gradT = gradT/dl;

        // Compute diffusion coefficient at the interface
        real kappa;
        if(haveThConductivity == UserDefFunction) {
          kappa = HALF_F*(kappaArr(k,j,i) +  kappaArr(k-koffset,j-joffset,i-ioffset));
        } else {
          kappa = kappaConstant;
        }

        // Add thermal diffusion to the flux
        Flux(ENG,k,j,i) += -kappa*gradT;

        // Compute total diffusion coefficient
        real locdmax = kappa * gamma_m1 /
                        (HALF_F * ( Vc(RHO,k,j,i) + Vc(RHO,k-koffset,j-joffset,i-ioffset)));
        dMax(k,j,i) = FMAX(dMax(k,j,i) , locdmax);
      });
  } else { //Braginskii
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
    real knorConstant = this->knor;
    real kparConstant = this->kpar;
    IdefixArray3D<real> knorArr = this->knorArr;
    IdefixArray3D<real> kparArr = this->kparArr;

    if(haveThConductivity == UserDefFunction && dir == IDIR) {
      if(bragThermalConductivityFunc) {
        idfx::pushRegion("UserDef::bragThermalDiffusivityFunction");
        bragThermalConductivityFunc(*this->hydro->data, t, kparArr, knorArr);
        idfx::popRegion();
      } else {
        IDEFIX_ERROR("No user-defined thermal diffusion function has been enrolled");
      }
    }

    idefix_for("BragThermalFlux",kbeg, kend, jbeg, jend, ibeg, iend,
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
          if(haveThConductivity == UserDefFunction) {
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

          if(haveThConductivity == UserDefFunction) {
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
//#                dTk = D_DZ_J_tmp(0,RHO)/dx3(k);
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
          if(haveThConductivity == UserDefFunction) {
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

        dT_mag = SQRT(dT_mag) + 0.01*SMALL_NUMBER;

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
        Bmag = SQRT(Bmag);
        Bmag = FMAX(1e-6*SMALL_NUMBER,Bmag);

        bgradT /= Bmag;

        bn = Bn/Bmag; /* -- unit vector component -- */
        q = kpar*bgradT*bn + knor*(dTn - bn*bgradT);
        //qi_mag = SQRT(  (kpar*kpar - knor*knor)*bgradT*bgradT
        //  + knor*knor*dT_mag*dT_mag);

        Flux(ENG, k, j, i) -= q;

        dMax(k,j,i) += FMAX(dMax(k,j,i),locdmax);
      });
  } // End of Braginskii

  idfx::popRegion();
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

#endif // HAVE_ENERGY
