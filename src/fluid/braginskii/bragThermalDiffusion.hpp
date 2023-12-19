// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_BRAGINSKII_BRAGTHERMALDIFFUSION_HPP_
#define FLUID_BRAGINSKII_BRAGTHERMALDIFFUSION_HPP_

#include <string>

#include "idefix.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "fluid_defs.hpp"
#include "eos.hpp"
#include "slopeLimiter.hpp"

// Forward class hydro declaration
template <typename Phys> class Fluid;

class DataBlock;

class BragThermalDiffusion {
 public:
  template <typename Phys>
  BragThermalDiffusion(Input &, Grid &, Fluid<Phys> *);  // Initialisation

  void ShowConfig(); // display configuration

  void AddBragDiffusiveFlux(int, const real, const IdefixArray4D<real> &);

  template<const PLMLimiter>
  void AddBragDiffusiveFluxLim(int, const real, const IdefixArray4D<real> &);

  // Enroll user-defined thermal conductivity
  void EnrollBragThermalDiffusivity(BragDiffusivityFunc);

  IdefixArray3D<real> heatSrc;  // Source terms of the thermal operator
  IdefixArray3D<real> knorArr;
  IdefixArray3D<real> kparArr;

  // pre-computed geometrical factors in non-cartesian geometry
  IdefixArray1D<real> one_dmu;

 private:
  DataBlock *data;

  // status of the module
  ParabolicModuleStatus &status;

  BragDiffusivityFunc diffusivityFunc;

  bool haveSlopeLimiter{false};

  // helper array
  IdefixArray4D<real> &Vc;
  IdefixArray4D<real> &Vs;
  IdefixArray3D<real> &dMax;

  // constant diffusion coefficient (when needed)
  real knor, kpar;

  // equation of state (required to get the heat capacity)
  EquationOfState *eos;

  PLMLimiter limiter{PLMLimiter::VanLeer};
};

#include "fluid.hpp"
#include "dataBlock.hpp"

template <typename Phys>
BragThermalDiffusion::BragThermalDiffusion(Input &input, Grid &grid, Fluid<Phys> *hydroin):
                            Vc(hydroin->Vc),
                            Vs(hydroin->Vs),
                            dMax(hydroin->dMax),
                            eos(hydroin->eos.get()),
                            data(hydroin->data),
                            status(hydroin->bragThermalDiffusionStatus) {
  idfx::pushRegion("BragThermalDiffusion::BragThermalDiffusion");

  if(input.CheckEntry("Hydro","bragTDiffusion")>=0) {
    if(input.Get<std::string>("Hydro","bragTDiffusion",1).compare("mc") == 0) {
      this->haveSlopeLimiter = true;
      limiter = PLMLimiter::McLim;
    } else if(input.Get<std::string>("Hydro","bragTDiffusion",1).compare("vanleer") == 0) {
      this->haveSlopeLimiter = true;
      limiter = PLMLimiter::VanLeer;
    } else if(input.Get<std::string>("Hydro","bragTDiffusion",1).compare("minmod") == 0) {
      IDEFIX_ERROR("The minmod slope limiter is not available because it has been "
                   "found to be too diffusive.");
    } else if(input.Get<std::string>("Hydro","bragTDiffusion",1).compare("nolimiter") == 0) {
      this->haveSlopeLimiter = false;
//      limiter = PLMLimiter::VanLeer;
    } else {
      IDEFIX_ERROR("Unknown braginskii thermal diffusion limiter in idefix.ini. "
                   "Can only be vanleer, mc or nolimiter.");
    }
    if(input.Get<std::string>("Hydro","bragTDiffusion",2).compare("constant") == 0) {
      this->kpar = input.Get<real>("Hydro","bragTDiffusion",3);
      this->knor = input.GetOrSet<real>("Hydro","bragTDiffusion",4,0.);
      this->status.status = Constant;
    } else if(input.Get<std::string>("Hydro","bragTDiffusion",2).compare("userdef") == 0) {
      this->status.status = UserDefFunction;
      this->kparArr = IdefixArray3D<real>("BragThermalDiffusionKparArray",data->np_tot[KDIR],
                                                               data->np_tot[JDIR],
                                                               data->np_tot[IDIR]);
      this->knorArr = IdefixArray3D<real>("BragThermalDiffusionKnorArray",data->np_tot[KDIR],
                                                               data->np_tot[JDIR],
                                                               data->np_tot[IDIR]);
    } else {
      IDEFIX_ERROR("Unknown braginskii thermal diffusion definition in idefix.ini. "
                   "Can only be constant or userdef.");
    }
  } else {
    IDEFIX_ERROR("I cannot create a BragThermalDiffusion object without bragTDiffusion defined"
                   "in the .ini file");
  }

  #ifndef MHD
    IDEFIX_ERROR("Braginskii Thermal diffusion requires MHD");
  #endif
  #ifdef ISOTHERMAL
    IDEFIX_ERROR("Braginskii Thermal diffusion is not compatible"
                 "with the ISOTHERMAL approximation");
  #endif

  idfx::popRegion();
}

//We now define spatial derivative macros for the temperature field.
//    The temperature is not a primitive variable
//    and therefore not available as such in the DataBlock.
//    It is rather defined as PRS/RHO.
//    Special spatial derivative macros are therefore needed and defined here
//    directly at the right cell interface according to the direciton of the flux.
#define D_DX_I_T(q)  (q(PRS,k,j,i)/q(RHO,k,j,i) - q(PRS,k,j,i - 1)/q(RHO,k,j,i - 1))
#define D_DY_J_T(q)  (q(PRS,k,j,i)/q(RHO,k,j,i) - q(PRS,k,j - 1,i)/q(RHO,k,j - 1,i))
#define D_DZ_K_T(q)  (q(PRS,k,j,i)/q(RHO,k,j,i) - q(PRS,k - 1,j,i)/q(RHO,k - 1,j,i))

#define SL_DX_T(q,k,j,i)  (q(PRS,k,j,i)/q(RHO,k,j,i)           \
                                        - q(PRS,k,j,i - 1)/q(RHO,k,j,i - 1))
#define SL_DY_T(q,k,j,i)  (q(PRS,k,j,i)/q(RHO,k,j,i)            \
                                        - q(PRS,k,j - 1,i)/q(RHO,k,j - 1,i))
#define SL_DZ_T(q,k,j,i)  (q(PRS,k,j,i)/q(RHO,k,j,i)            \
                                        - q(PRS,k - 1,j,i)/q(RHO,k - 1,j,i))

#define D_DY_I_T(q)  (  0.25*(q(PRS,k,j + 1,i    ) / q(RHO,k,j + 1,i)          \
                                  + q(PRS,k,j + 1,i - 1) / q(RHO,k,j + 1,i - 1))     \
                            - 0.25*(q(PRS,k,j - 1,i)     / q(RHO,k,j - 1,i)          \
                                  + q(PRS,k,j - 1,i - 1) / q(RHO,k,j - 1,i - 1)))

#define D_DZ_I_T(q)  (  0.25*(q(PRS,k + 1,j,i)     / q(RHO,k + 1,j,i)       \
                                  + q(PRS,k + 1,j,i - 1) / q(RHO,k + 1,j,i - 1))  \
                            - 0.25*(q(PRS,k - 1,j,i)     / q(RHO,k - 1,j,i)       \
                                  + q(PRS,k - 1,j,i - 1) / q(RHO,k - 1,j,i - 1)))

#define D_DX_J_T(q)  (  0.25*(q(PRS,k,j,i + 1)     / q(RHO,k,j,i + 1)      \
                                  + q(PRS,k,j - 1,i + 1) / q(RHO,k,j - 1,i + 1)) \
                            - 0.25*(q(PRS,k,j,i - 1)     / q(RHO,k,j,i - 1)      \
                                  + q(PRS,k,j - 1,i - 1) / q(RHO,k,j - 1,i - 1)))

#define D_DZ_J_T(q)  (  0.25*(q(PRS,k + 1,j,i)     / q(RHO,k + 1,j,i)       \
                                  + q(PRS,k + 1,j - 1,i) / q(RHO,k + 1,j - 1,i))  \
                            - 0.25*(q(PRS,k - 1,j,i)     / q(RHO,k - 1,j,i)       \
                                  + q(PRS,k - 1,j - 1,i) / q(RHO,k - 1,j - 1,i)))

#define D_DX_K_T(q)  (  0.25*(q(PRS,k,j,i + 1)     / q(RHO,k,j,i + 1)      \
                                  + q(PRS,k - 1,j,i + 1) / q(RHO,k - 1,j,i + 1)) \
                            - 0.25*(q(PRS,k,j,i - 1)     / q(RHO,k,j,i - 1)      \
                                  + q(PRS,k - 1,j,i - 1) / q(RHO,k - 1,j,i - 1)))

#define D_DY_K_T(q)  (  0.25*(q(PRS,k,j + 1,i)     / q(RHO,k,j + 1,i)       \
                                  + q(PRS,k - 1,j + 1,i) / q(RHO,k - 1,j + 1,i))  \
                            - 0.25*(q(PRS,k,j - 1,i)     / q(RHO,k,j - 1,i)       \
                                  + q(PRS,k - 1,j - 1,i) / q(RHO,k - 1,j - 1,i))) \

//We now define spatial average macros for the magnetic field.
//    The magnetic field appears in the expression of the Braginskii heat flux.
//    It is therefore needed at the right cell interface according to the direction of the flux.
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


// (this avoids an extra array)
template <PLMLimiter limTemplate>
void BragThermalDiffusion::AddBragDiffusiveFluxLim(int dir, const real t,
                                                const IdefixArray4D<real> &Flux) {
  idfx::pushRegion("BragThermalDiffusion::AddBragDiffusiveFluxLim");

  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Vs = this->Vs;
  IdefixArray3D<real> dMax = this->dMax;
  EquationOfState eos = *(this->eos);

  HydroModuleStatus haveThermalDiffusion = this->status.status;
  bool haveSlopeLimiter = this->haveSlopeLimiter;
  using SL = SlopeLimiter<limTemplate>;

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

  #if GEOMETRY == SPHERICAL
    IdefixArray1D<real> rt   = this->data->rt;
    IdefixArray1D<real> dmu  = this->data->dmu;
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
      real q;
      [[maybe_unused]] real Bi, Bj, Bk, Bn;
      Bi = Bj = Bk = Bn = ZERO_F;

      [[maybe_unused]] real dTi, dTj, dTk, dTn;
      dTi = dTj = dTk = dTn = ZERO_F;

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

        EXPAND( Bi = BX_I; ,
                Bj = BY_I; ,
                Bk = BZ_I; )
        Bn = BX_I;

        #if GEOMETRY == CARTESIAN
          dTi = D_DX_I_T(Vc)/dx1(i);
          #if DIMENSIONS >= 2
            if (haveSlopeLimiter) {
              dTj = SL::PLMLim(SL_DY_T(Vc,k,j,i-1)/dx2(j),
                                 SL_DY_T(Vc,k,j+1,i-1)/dx2(j+1));
              dTj = SL::PLMLim(dTj,
                                 SL::PLMLim(SL_DY_T(Vc,k,j,i)/dx2(j),
                                              SL_DY_T(Vc,k,j+1,i)/dx2(j+1)));
            } else {
              dTj = D_DY_I_T(Vc)/dx2(j);
            }
            #if DIMENSIONS == 3
              if (haveSlopeLimiter) {
                dTk = SL::PLMLim(SL_DZ_T(Vc,k,j,i-1)/dx3(k),
                                   SL_DZ_T(Vc,k+1,j,i-1)/dx3(k+1));
                dTk = SL::PLMLim(dTk,
                                   SL::PLMLim(SL_DZ_T(Vc,k,j,i)/dx3(k),
                                                SL_DZ_T(Vc,k+1,j,i)/dx3(k+1)));
              } else {
                dTk = D_DZ_I_T(Vc)/dx3(k);
              }
            #endif
          #endif
        #elif GEOMETRY == CYLINDRICAL
          dTi = D_DX_I_T(Vc)/dx1(i);
          #if DIMENSIONS >= 2
            if (haveSlopeLimiter) {
              dTj = SL::PLMLim(SL_DY_T(Vc,k,j,i-1)/dx2(j),
                                 SL_DY_T(Vc,k,j+1,i-1)/dx2(j+1));
              dTj = SL::PLMLim(dTj,
                                 SL::PLMLim(SL_DY_T(Vc,k,j,i)/dx2(j),
                                              SL_DY_T(Vc,k,j+1,i)/dx2(j+1)));
            } else {
              dTj = D_DY_I_T(Vc)/dx2(j);
            }
          #endif
          // No cylindrical geometry in 3D!
        #elif GEOMETRY == POLAR
          dTi = D_DX_I_T(Vc)/dx1(i);
          #if DIMENSIONS >= 2
            if (haveSlopeLimiter) {
              dTj = 1./x1(i-1)*SL::PLMLim(SL_DY_T(Vc,k,j,i-1)/dx2(j),
                                            SL_DY_T(Vc,k,j+1,i-1)/dx2(j+1));
              dTj = SL::PLMLim(dTj,
                                 1./x1(i)*SL::PLMLim(SL_DY_T(Vc,k,j,i)/dx2(j),
                                                       SL_DY_T(Vc,k,j+1,i)/dx2(j+1)));
            } else {
              dTj = 1./x1l(i)*D_DY_I_T(Vc)/dx2(j); // 1/r dTj/dxj
            }
            #if DIMENSIONS == 3
              if (haveSlopeLimiter) {
                dTk = SL::PLMLim(SL_DZ_T(Vc,k,j,i-1)/dx3(k),
                                   SL_DZ_T(Vc,k+1,j,i-1)/dx3(k+1));
                dTk = SL::PLMLim(dTk,
                                   SL::PLMLim(SL_DZ_T(Vc,k,j,i)/dx3(k),
                                                SL_DZ_T(Vc,k+1,j,i)/dx3(k+1)));
              } else {
                dTk = D_DZ_I_T(Vc)/dx3(k);
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

          dTi = D_DX_I_T(Vc)/dx1(i);
          #if DIMENSIONS >= 2
            if (haveSlopeLimiter) {
              dTj = 1./x1(i-1)*SL::PLMLim(SL_DY_T(Vc,k,j,i-1)/dx2(j),
                                            SL_DY_T(Vc,k,j+1,i-1)/dx2(j+1));
              dTj = SL::PLMLim(dTj,
                    1./x1(i)*SL::PLMLim(SL_DY_T(Vc,k,j,i)/dx2(j),
                                          SL_DY_T(Vc,k,j+1,i)/dx2(j+1)));
            } else {
              dTj = 1./x1l(i)*D_DY_I_T(Vc)/dx2(j); // 1/r dTj/dxj
            }
            #if DIMENSIONS == 3
              if (haveSlopeLimiter) {
                dTk = s_1/x1(i-1)*SL::PLMLim(SL_DZ_T(Vc,k,j,i-1)/dx3(k),
                                               SL_DZ_T(Vc,k+1,j,i-1)/dx3(k+1));
                dTk = SL::PLMLim(dTk,
                                   s_1/x1(i)*SL::PLMLim(SL_DZ_T(Vc,k,j,i)/dx3(k),
                                                          SL_DZ_T(Vc,k+1,j,i)/dx3(k+1)));
              } else {
                dTk = s_1/x1l(i)*D_DZ_I_T(Vc)/dx3(k);
              }
            #endif
          #endif
        #endif // GEOMETRY

        real gamma_m1 = eos.GetGamma(Vc(PRS,k,j,i),Vc(RHO,k,j,i)) - ONE_F;

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

        EXPAND( Bi = BX_J; ,
                Bj = BY_J; ,
                Bk = BZ_J; )
        Bn = BY_J;

        #if GEOMETRY == CARTESIAN
          if (haveSlopeLimiter) {
            dTi = SL::PLMLim(SL_DX_T(Vc,k,j-1,i)/dx1(i),
                              SL_DX_T(Vc,k,j-1,i+1)/dx1(i+1));
            dTi = SL::PLMLim(dTi,
                               SL::PLMLim(SL_DX_T(Vc,k,j,i)/dx1(i),
                                            SL_DX_T(Vc,k,j,i+1)/dx1(i+1)));
          } else {
            dTi = D_DX_J_T(Vc)/dx1(i);
          }
          dTj = D_DY_J_T(Vc)/dx2(j);
          #if DIMENSIONS == 3
            if (haveSlopeLimiter) {
              dTk = SL::PLMLim(SL_DX_T(Vc,k,j-1,i)/dx3(k),
                                 SL_DX_T(Vc,k+1,j-1,i)/dx3(k+1));
              dTk = SL::PLMLim(dTk,
                                 SL::PLMLim(SL_DX_T(Vc,k,j,i)/dx3(k),
                                              SL_DX_T(Vc,k+1,j,i)/dx3(k+1)));
            } else {
              dTk = D_DZ_J_T(Vc)/dx3(k);
            }
          #endif
        #elif GEOMETRY == CYLINDRICAL
          if (haveSlopeLimiter) {
            dTi = SL::PLMLim(SL_DX_T(Vc,k,j-1,i)/dx1(i),
                               SL_DX_T(Vc,k,j-1,i+1)/dx1(i+1));
            dTi = SL::PLMLim(dTi,
                               SL::PLMLim(SL_DX_T(Vc,k,j,i)/dx1(i),
                                            SL_DX_T(Vc,k,j,i+1)/dx1(i+1)));
          } else {
            dTi = D_DX_J_T(Vc)/dx1(i);
          }
          dTj = D_DY_J_T(Vc)/dx2(j);
          // No cylindrical geometry in 3D!
        #elif GEOMETRY == POLAR
          //gradT = ... + 1/r dT/dtheta + ...
          if (haveSlopeLimiter) {
            dTi = SL::PLMLim(SL_DX_T(Vc,k,j-1,i)/dx1(i),
                               SL_DX_T(Vc,k,j-1,i+1)/dx1(i+1));
            dTi = SL::PLMLim(dTi,
                               SL::PLMLim(SL_DX_T(Vc,k,j,i)/dx1(i),
                                            SL_DX_T(Vc,k,j,i+1)/dx1(i+1)));
          } else {
            dTi = D_DX_J_T(Vc)/dx1(i);
          }
          dTj = 1./x1(i)*D_DY_J_T(Vc)/dx2(j);
          #if DIMENSIONS == 3
            if (haveSlopeLimiter) {
              dTk = SL::PLMLim(SL_DX_T(Vc,k,j-1,i)/dx3(k),
                                 SL_DX_T(Vc,k+1,j-1,i)/dx3(k+1));
              dTk = SL::PLMLim(dTk,
                                 SL::PLMLim(SL_DX_T(Vc,k,j,i)/dx3(k),
                                              SL_DX_T(Vc,k+1,j,i)/dx3(k+1)));
            } else {
              dTk = D_DZ_J_T(Vc)/dx3(k);
            }
          #endif
        #elif GEOMETRY == SPHERICAL
          if (haveSlopeLimiter) {
            dTi = SL::PLMLim(SL_DX_T(Vc,k,j-1,i)/dx1(i),
                               SL_DX_T(Vc,k,j-1,i+1)/dx1(i+1));
            dTi = SL::PLMLim(dTi,
                               SL::PLMLim(SL_DX_T(Vc,k,j,i)/dx1(i),
                                            SL_DX_T(Vc,k,j,i+1)/dx1(i+1)));
          } else {
            dTi = D_DX_J_T(Vc)/dx1(i);
          }
          dTj = 1./x1(i)*D_DY_J_T(Vc)/dx2(j);
          #if DIMENSIONS == 3
            if (haveSlopeLimiter) {
              real s_1 = sinx2(j-1);

              // Trick to ensure that the axis does not lead to Nans
              if(FABS(s_1) < SMALL_NUMBER) {
                s_1 = ZERO_F;
              } else {
                s_1 = ONE_F/s_1;
              }
              dTk = s_1/x1(i)*SL::PLMLim(SL_DX_T(Vc,k,j-1,i)/dx3(k),
                                           SL_DX_T(Vc,k+1,j-1,i)/dx3(k+1));

              s_1 = sinx2(j);

              // Trick to ensure that the axis does not lead to Nans
              if(FABS(s_1) < SMALL_NUMBER) {
                s_1 = ZERO_F;
              } else {
                s_1 = ONE_F/s_1;
              }
              dTk = s_1/x1(i)*SL::PLMLim(dTk,
                                           SL::PLMLim(SL_DX_T(Vc,k,j,i)/dx3(k),
                                                        SL_DX_T(Vc,k+1,j,i)/dx3(k+1)));
            } else {
              real s_1 = sinx2m(j);

              // Trick to ensure that the axis does not lead to Nans
              if(FABS(s_1) < SMALL_NUMBER) {
                s_1 = ZERO_F;
              } else {
                s_1 = ONE_F/s_1;
              }

              dTk = s_1/x1(i)*D_DZ_J_T(Vc)/dx3(k);
            }
          #endif
        #endif // GEOMETRY

        real gamma_m1 = eos.GetGamma(Vc(PRS,k,j,i),Vc(RHO,k,j,i)) - ONE_F;

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
            dTi = SL::PLMLim(SL_DX_T(Vc,k-1,j,i)/dx1(i),
                               SL_DX_T(Vc,k-1,j,i+1)/dx1(i+1));
            dTi = SL::PLMLim(dTi,
                               SL::PLMLim(SL_DX_T(Vc,k,j,i)/dx1(i),
                                            SL_DX_T(Vc,k,j,i+1)/dx1(i+1)));
            dTj = SL::PLMLim(SL_DY_T(Vc,k-1,j,i)/dx2(j),
                               SL_DY_T(Vc,k-1,j+1,i)/dx2(j+1));
            dTj = SL::PLMLim(dTj,
                               SL::PLMLim(SL_DY_T(Vc,k,j,i)/dx2(j),
                                            SL_DY_T(Vc,k,j+1,i)/dx2(j+1)));
          } else {
            dTi = D_DX_K_T(Vc)/dx1(i);
            dTj = D_DY_K_T(Vc)/dx2(j);
          }
          dTk = D_DZ_K_T(Vc)/dx3(k);
        // No cylindrical geometry in 3D!
        #elif GEOMETRY == POLAR
          if (haveSlopeLimiter) {
            dTi = SL::PLMLim(SL_DX_T(Vc,k-1,j,i)/dx1(i),
                               SL_DX_T(Vc,k-1,j,i+1)/dx1(i+1));
            dTi = SL::PLMLim(dTi,
                               SL::PLMLim(SL_DX_T(Vc,k,j,i)/dx1(i),
                                            SL_DX_T(Vc,k,j,i+1)/dx1(i+1)));
            dTj = SL::PLMLim(SL_DY_T(Vc,k-1,j,i)/dx2(j),
                               SL_DY_T(Vc,k-1,j+1,i)/dx2(j+1));
            dTj = SL::PLMLim(dTj,
                               SL::PLMLim(SL_DY_T(Vc,k,j,i)/dx2(j),
                                            SL_DY_T(Vc,k,j+1,i)/dx2(j+1)));
            dTj *= 1./x1(i);
          } else {
            dTi = D_DX_K_T(Vc)/dx1(i);
            dTj = 1./x1(i)*D_DY_K_T(Vc)/dx2(j); // 1/r dTj/dxj
          }
          dTk = D_DZ_K_T(Vc)/dx3(k);
        #elif GEOMETRY == SPHERICAL
            real s_1 = sinx2(j);

          // Trick to ensure that the axis does not lead to Nans
          if(FABS(s_1) < SMALL_NUMBER) {
            s_1 = ZERO_F;
          } else {
            s_1 = ONE_F/s_1;
          }

          if (haveSlopeLimiter) {
            dTi = SL::PLMLim(SL_DX_T(Vc,k-1,j,i)/dx1(i),
                               SL_DX_T(Vc,k-1,j,i+1)/dx1(i+1));
            dTi = SL::PLMLim(dTi,
                               SL::PLMLim(SL_DX_T(Vc,k,j,i)/dx1(i),
                                            SL_DX_T(Vc,k,j,i+1)/dx1(i+1)));
            dTj = SL::PLMLim(SL_DY_T(Vc,k-1,j,i)/dx2(j),
                               SL_DY_T(Vc,k-1,j+1,i)/dx2(j+1));
            dTj = SL::PLMLim(dTj,
                               SL::PLMLim(SL_DY_T(Vc,k,j,i)/dx2(j),
                                            SL_DY_T(Vc,k,j+1,i)/dx2(j+1)));
            dTj *= 1./x1(i);
          } else {
            dTi = D_DX_K_T(Vc)/dx1(i);
            dTj = 1./x1(i)*D_DY_K_T(Vc)/dx2(j); // 1/r dTj/dxj
          }
          dTk = s_1/x1(i)*D_DZ_K_T(Vc)/dx3(k);
          //gradT = ... + ... + 1/(r*sin(theta)) dTphi/dphi
        #endif // GEOMETRY

        real gamma_m1 = eos.GetGamma(Vc(PRS,k,j,i),Vc(RHO,k,j,i)) - ONE_F;

        locdmax = FMAX(kpar,knor)/(0.5*(Vc(RHO,k,j,i) + Vc(RHO,k-1,j,i)))*gamma_m1;

        dTn = dTk;
      }
      // From here, gradients and normal have been computed, so we just need to get the fluxes

      bgradT = EXPAND( Bi*dTi , + Bj*dTj, +Bk*dTk);
      Bmag = EXPAND( Bi*Bi , + Bj*Bj, + Bk*Bk);
      Bmag = sqrt(Bmag);
      Bmag = FMAX(1e-6*SMALL_NUMBER,Bmag);

      bgradT /= Bmag;

      bn = Bn/Bmag; /* -- unit vector component -- */
      q = kpar*bgradT*bn + knor*(dTn - bn*bgradT);

      Flux(ENG, k, j, i) -= q;

      dMax(k,j,i) = FMAX(dMax(k,j,i),locdmax);
    });
  idfx::popRegion();
}

#undef D_DX_I_T
#undef D_DY_J_T
#undef D_DZ_K_T
#undef SL_DX_T
#undef SL_DY_T
#undef SL_DZ_T
#undef D_DY_I_T
#undef D_DZ_I_T
#undef D_DX_J_T
#undef D_DZ_J_T
#undef D_DX_K_T
#undef D_DY_K_T
#undef BX_I
#undef BX_J
#undef BX_K
#undef BY_I
#undef BY_J
#undef BY_K
#undef BZ_I
#undef BZ_J
#undef BZ_K
#endif // FLUID_BRAGINSKII_BRAGTHERMALDIFFUSION_HPP_
