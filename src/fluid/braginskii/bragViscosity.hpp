// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_BRAGINSKII_BRAGVISCOSITY_HPP_
#define FLUID_BRAGINSKII_BRAGVISCOSITY_HPP_

#include <string>

#include "idefix.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "fluid_defs.hpp"
#include "slopeLimiter.hpp"

// Forward class hydro declaration
template <typename Phys> class Fluid;
class DataBlock;

class BragViscosity {
 public:
  template <typename Phys>
  BragViscosity(Input &, Grid &, Fluid<Phys> *);
  void ShowConfig();                    // print configuration
  void AddBragViscousFlux(int, const real, const IdefixArray4D<real> &);

  template <const PLMLimiter>
  void AddBragViscousFluxLim(int, const real, const IdefixArray4D<real> &);

  // Enroll user-defined viscous diffusivity
  void EnrollBragViscousDiffusivity(DiffusivityFunc);

  // Function for internal use (but public to allow for Cuda lambda capture)
  void InitArrays();

  IdefixArray4D<real> bragViscSrc;  // Source terms of the viscous operator
  IdefixArray3D<real> etaBragArr;


  // pre-computed geometrical factors in non-cartesian geometry
  IdefixArray1D<real> one_dmu;

 private:
  DataBlock* data;

  // Viscosity status
  ParabolicModuleStatus &status;

  DiffusivityFunc bragViscousDiffusivityFunc;

  bool haveSlopeLimiter{false};

  IdefixArray4D<real> &Vc;
  IdefixArray4D<real> &Vs;
  IdefixArray3D<real> &dMax;

  // constant diffusion coefficient (when needed)
  real etaBrag;

  PLMLimiter limiter{PLMLimiter::VanLeer};
};

#include "fluid.hpp"

template<typename Phys>
BragViscosity::BragViscosity(Input &input, Grid &grid, Fluid<Phys> *hydroin):
                      Vc(hydroin->Vc),
                      Vs(hydroin->Vs),
                      dMax(hydroin->dMax),
                      status(hydroin->bragViscosityStatus) {
  idfx::pushRegion("BragViscosity::BragViscosity");
  // Save the parent hydro object
  this->data = hydroin->data;

  if(input.CheckEntry("Hydro","bragViscosity")>=0) {
    if(input.Get<std::string>("Hydro","bragViscosity",1).compare("vanleer") == 0) {
      this->haveSlopeLimiter = true;
      limiter = PLMLimiter::VanLeer;
    } else if(input.Get<std::string>("Hydro","bragViscosity",1).compare("mc") == 0) {
      this->haveSlopeLimiter = true;
      limiter = PLMLimiter::McLim;
    } else if (input.Get<std::string>("Hydro","bragViscosity",1).compare("nolimiter") == 0) {
      this->haveSlopeLimiter = false;
//      limiter = PLMLimiter::VanLeer;
    } else {
      IDEFIX_ERROR("Unknown braginskii viscosity limiter in idefix.ini. "
                   "Can only be vanleer, mc or nolimiter.");
    }
    if(input.Get<std::string>("Hydro","bragViscosity",2).compare("constant") == 0) {
        this->etaBrag = input.Get<real>("Hydro","bragViscosity",3);
        this->status.status = Constant;
    } else if(input.Get<std::string>("Hydro","bragViscosity",2).compare("userdef") == 0) {
        this->status.status = UserDefFunction;
        this->etaBragArr = IdefixArray3D<real>("BragViscosityEtaArray",data->np_tot[KDIR],
                                                                 data->np_tot[JDIR],
                                                                 data->np_tot[IDIR]);
    } else {
      IDEFIX_ERROR("Unknown braginskii viscosity definition in idefix.ini. "
                   "Can only be constant or userdef.");
    }
  } else {
    IDEFIX_ERROR("I cannot create a BragViscosity object without viscosity defined"
                   "in the .ini file");
  }

  InitArrays();

  idfx::popRegion();
}


//We now define spatial derivative macros for the velocity field.
//    They are directly computed at the right cell interface according to the direciton of the flux.
#define D_DX_I(q,n)  (q(n,k,j,i) - q(n,k,j,i - 1))
#define D_DY_J(q,n)  (q(n,k,j,i) - q(n,k,j - 1,i))
#define D_DZ_K(q,n)  (q(n,k,j,i) - q(n,k - 1,j,i))

#define SL_DX(q,n,k,j,i)  (q(n,k,j,i) - q(n,k,j,i - 1))
#define SL_DY(q,n,k,j,i)  (q(n,k,j,i) - q(n,k,j - 1,i))
#define SL_DZ(q,n,k,j,i)  (q(n,k,j,i) - q(n,k - 1,j,i))

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


// This function computes the Braginskii viscous flux and stores it in hydro->fluxRiemann
// (this avoids an extra array)
// Associated source terms, present in non-cartesian geometry are also computed
// and stored in this->bragViscSrc for later use (in calcRhs).
template <PLMLimiter limTemplate>
void BragViscosity::AddBragViscousFluxLim(int dir, const real t, const IdefixArray4D<real> &Flux) {
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
  bool haveSlopeLimiter = this->haveSlopeLimiter;

  using SL = SlopeLimiter<limTemplate>;

  // Braginskii Viscosity
  // Compute viscosity if needed
  if(haveViscosity == UserDefFunction && dir == IDIR) {
    if(bragViscousDiffusivityFunc) {
      bragViscousDiffusivityFunc(*this->data, t, etaBragArr);
    } else {
      IDEFIX_ERROR("No user-defined Braginskii viscosity function has been enrolled");
    }
  }

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

      [[maybe_unused]] real tau_xx, tau_xy, tau_xz;
      [[maybe_unused]] real tau_yy, tau_yz;
      [[maybe_unused]] real tau_zz;

      [[maybe_unused]] real dVxi, dVxj, dVxk;
      [[maybe_unused]] real dVyi, dVyj, dVyk;
      [[maybe_unused]] real dVzi, dVzj, dVzk;

      real divV;

      bi = bj = bk = ZERO_F;

      tau_xx = tau_xy = tau_xz = ZERO_F;
      tau_yy = tau_yz = ZERO_F;
      tau_zz = ZERO_F;

      dVxi = dVxj = dVxk = ZERO_F;
      dVyi = dVyj = dVyk = ZERO_F;
      dVzi = dVzj = dVzk = ZERO_F;

      real etaBrag;
      [[maybe_unused]] real etaBragC;


      //Cell-centered values to compute source terms
      #if GEOMETRY != CARTESIAN
        real biC, bjC, bkC;
        real BmagC;
        real Pnor_parC;
        real bbgradVC;

        [[maybe_unused]] real tau_xxC, tau_xyC, tau_xzC;
        [[maybe_unused]] real tau_yyC, tau_yzC;
        [[maybe_unused]] real tau_zzC;

        [[maybe_unused]] real dVxiC, dVxjC, dVxkC;
        [[maybe_unused]] real dVyiC, dVyjC, dVykC;
        [[maybe_unused]] real dVziC, dVzjC, dVzkC;

        real divVC;

        biC = bjC = bkC = ZERO_F;
        tau_xxC = tau_xyC = tau_xzC = ZERO_F;
        tau_yyC = tau_yzC = ZERO_F;
        tau_zzC = ZERO_F;

        dVxiC = dVxjC = dVxkC = ZERO_F;
        dVyiC = dVyjC = dVykC = ZERO_F;
        dVziC = dVzjC = dVzkC = ZERO_F;

        //Compute values at the center of the cells, no slope limiter are used at this stage
        EXPAND(  dVxiC = HALF_F*(Vc(VX1,k,j,i + 1) - Vc(VX1,k,j,i - 1))/dx1(i); ,
                 dVyiC = HALF_F*(Vc(VX2,k,j,i + 1) - Vc(VX2,k,j,i - 1))/dx1(i); ,
                 dVziC = HALF_F*(Vc(VX3,k,j,i + 1) - Vc(VX3,k,j,i - 1))/dx1(i); )
        #if DIMENSIONS >= 2
          EXPAND(  dVxjC = HALF_F*(Vc(VX1,k,j + 1,i) - Vc(VX1,k,j - 1,i))/dx2(j); ,
                   dVyjC = HALF_F*(Vc(VX2,k,j + 1,i) - Vc(VX2,k,j - 1,i))/dx2(j); ,
                   dVzjC = HALF_F*(Vc(VX3,k,j + 1,i) - Vc(VX3,k,j - 1,i))/dx2(j); )
          #if DIMENSIONS == 3
            EXPAND (  dVxkC = HALF_F*(Vc(VX1,k + 1,j,i) - Vc(VX1,k - 1,j,i))/dx3(k); ,
                      dVykC = HALF_F*(Vc(VX2,k + 1,j,i) - Vc(VX2,k - 1,j,i))/dx3(k); ,
                      dVzkC = HALF_F*(Vc(VX3,k + 1,j,i) - Vc(VX3,k - 1,j,i))/dx3(k); )
          #endif
        #endif

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
                     biC*biC*dVxiC + biC*bjC*dVyiC + biC*bkC*(dVziC  - Vc(VX3,k,j,i)/x1(i)) ,
                   + bjC*biC*dVxjC + bjC*bjC*dVyjC + bjC*bkC*dVzjC + bkC*bkC*Vc(VX1,k,j,i)/x1(i) ,
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
                                                                                + bkC*bkC*dVzkC
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
            etaBrag = 2.*(etaBragArr(k,j,i-1)*etaBragArr(k,j,i)) /
                         (etaBragArr(k,j,i-1)+etaBragArr(k,j,i));
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
            dVxj = SL::PLMLim(SL_DY(Vc,VX1,k,j + 1,i)/dx2(j+1),
                                SL_DY(Vc,VX1,k,j + 1,i - 1)/dx2(j+1));
            dVxj = SL::PLMLim(dVxj,
                                SL::PLMLim(SL_DY(Vc,VX1,k,j,i)/dx2(j),
                                             SL_DY(Vc,VX1,k,j,i - 1)/dx2(j)));
            dVyj = SL::PLMLim(SL_DY(Vc,VX2,k,j + 1,i)/dx2(j+1),
                                SL_DY(Vc,VX2,k,j + 1,i - 1)/dx2(j+1));
            dVyj = SL::PLMLim(dVyj,
                                SL::PLMLim(SL_DY(Vc,VX2,k,j,i)/dx2(j),
                                             SL_DY(Vc,VX2,k,j,i - 1)/dx2(j)));
            #if DIMENSIONS == 3
              dVzj = SL::PLMLim(SL_DY(Vc,VX3,k,j + 1,i)/dx2(j+1),
                                  SL_DY(Vc,VX3,k,j + 1,i - 1)/dx2(j+1));
              dVzj = SL::PLMLim(dVzj,
                                  SL::PLMLim(SL_DY(Vc,VX3,k,j,i)/dx2(j),
                                               SL_DY(Vc,VX3,k,j,i - 1)/dx2(j)));
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
              dVxk = SL::PLMLim(SL_DZ(Vc,VX1,k + 1,j,i)/dx3(k+1),
                                  SL_DZ(Vc,VX1,k + 1,j,i - 1)/dx3(k+1));
              dVxk = SL::PLMLim(dVxk,
                                  SL::PLMLim(SL_DZ(Vc,VX1,k,j,i)/dx3(k),
                                               SL_DZ(Vc,VX1,k,j,i - 1)/dx3(k)));
              dVyk = SL::PLMLim(SL_DZ(Vc,VX2,k + 1,j,i)/dx3(k+1),
                                  SL_DZ(Vc,VX2,k + 1,j,i - 1)/dx3(k+1));
              dVyk = SL::PLMLim(dVyk,
                                  SL::PLMLim(SL_DZ(Vc,VX2,k,j,i)/dx3(k),
                                               SL_DZ(Vc,VX2,k,j,i - 1)/dx3(k)));
              dVzk = SL::PLMLim(SL_DZ(Vc,VX3,k + 1,j,i)/dx3(k+1),
                                  SL_DZ(Vc,VX3,k + 1,j,i - 1)/dx3(k+1));
              dVzk = SL::PLMLim(dVzk,
                                  SL::PLMLim(SL_DZ(Vc,VX3,k,j,i)/dx3(k),
                                               SL_DZ(Vc,VX3,k,j,i - 1)/dx3(k)));
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

          divV = D_EXPAND(dVxi, + dVyj, + dVzk);
        // No source term in cartesian geometry
        #else
          [[maybe_unused]] real vx1i, vx2i, vx3i;
          vx1i = vx2i = vx3i = ZERO_F;

          EXPAND(  vx1i = 0.5*(Vc(VX1,k,j,i-1)+Vc(VX1,k,j,i)); ,
                   vx2i = 0.5*(Vc(VX2,k,j,i-1)+Vc(VX2,k,j,i)); ,
                   vx3i = 0.5*(Vc(VX3,k,j,i-1)+Vc(VX3,k,j,i)); )

          Pnor_parC = 3.0*etaBragC*(bbgradVC - divVC/3.0);
        #endif // GEOMETRY != CARTESIAN

        #if GEOMETRY == CYLINDRICAL
          bbgradV = D_EXPAND(
                       bi*bi*dVxi + bi*bj*dVyi + bi*bk*(dVzi  - vx3i/x1l(i)),
                     + bj*bi*dVxj + bj*bj*dVyj + bj*bk*dVzj + bk*bk*vx1i/x1l(i),
                       );

          divV = D_EXPAND(dVxi + vx1i/x1l(i) ,
                        + dVyj ,
                             );
          // No cylindrical geometry in 3D!

          //cell-centered values for source terms
          tau_zzC = Pnor_parC*(bkC*bkC - 1./3.);

          EXPAND( bragViscSrc(IDIR,k,j,i) = -tau_zzC/x1(i);  ,
                  bragViscSrc(JDIR,k,j,i) = ZERO_F;          ,
                  bragViscSrc(KDIR,k,j,i) = ZERO_F;          )

        #elif GEOMETRY == POLAR
          bbgradV = EXPAND(
                      bi*bi*dVxi                         + bj*bi*1./x1l(i)*dVxj
                                                                                  + bk*bi*dVxk,
                    + bi*bj*(dVyi - vx2i/x1l(i)) + bj*bj*(1./x1l(i)*dVyj + vx1i/x1l(i))
                                                                                  + bk*bj*dVyk,
                    + bi*bk*dVzi                         + bj*bk*1./x1l(i)*dVzj
                                                                                  + bk*bk*dVzk
                    );

          divV = D_EXPAND(dVxi + vx1i/x1l(i) ,
                           + 1./x1l(i)*dVyj ,
                           + dVzk);

          //cell-centered values for source terms
          tau_yyC = Pnor_parC*(bjC*bjC - 1./3.);

          EXPAND( bragViscSrc(IDIR,k,j,i) = -tau_yyC/x1(i);  ,
                  bragViscSrc(JDIR,k,j,i) = ZERO_F;          ,
                  bragViscSrc(KDIR,k,j,i) = ZERO_F;          )

        #elif GEOMETRY == SPHERICAL
          // NOLINT
          bbgradV = EXPAND(
      bi*bi*dVxi                 + bj*bi*1./x1l(i)*dVxj                 + bk*bi*s_1/x1l(i)*dVxk,
    + bi*bj*(dVyi - vx2i/x1l(i)) + bj*bj*(1./x1l(i)*dVyj + vx1i/x1l(i)) + bk*bj*s_1/x1l(i)*dVyk,
    + bi*bk*(dVzi - vx3i/x1l(i)) + bj*bk*(1./x1l(i)*dVzj - tan_1/x1l(i)*vx3i)
                                      + bk*bk*(s_1/x1l(i)*dVzk + vx1i/x1l(i) + tan_1/x1l(i)*vx2i));

          divV = D_EXPAND(2.0*vx1i/x1l(i) + dVxi,
                          + dVyj/x1l(i) + tan_1*vx2i/x1l(i),
                          + dVzk/x1l(i)*s_1 );

          //cell-centered values for source terms
          tau_yyC = Pnor_parC*(bjC*bjC - 1./3.);
          tau_zzC = Pnor_parC*(bkC*bkC - 1./3.);

          EXPAND( bragViscSrc(IDIR,k,j,i) = -(tau_yyC + tau_zzC)/x1(i);  ,
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
            etaBrag = 2.*(etaBragArr(k,j-1,i)*etaBragArr(k,j,i)) /
                         (etaBragArr(k,j-1,i)+etaBragArr(k,j,i));
          } else {
          etaBrag = HALF_F*(etaBragArr(k,j - 1,i)+etaBragArr(k,j,i));
          }
        } else {
          etaBragC = etaBrag = etaBragConstant;
        }

        if (haveSlopeLimiter) {
          dVxi = SL::PLMLim(SL_DX(Vc,VX1,k,j,i + 1)/dx1(i+1),
                              SL_DX(Vc,VX1,k,j - 1,i + 1)/dx1(i+1));
          dVxi = SL::PLMLim(dVxi,
                              SL::PLMLim(SL_DX(Vc,VX1,k,j,i)/dx1(i),
                                           SL_DX(Vc,VX1,k,j - 1,i)/dx1(i)));
          #if DIMENSIONS >= 2
            dVyi = SL::PLMLim(SL_DX(Vc,VX2,k,j,i + 1)/dx1(i+1),
                                SL_DX(Vc,VX2,k,j - 1,i + 1)/dx1(i+1));
            dVyi = SL::PLMLim(dVyi,
                                SL::PLMLim(SL_DX(Vc,VX2,k,j,i)/dx1(i),
                                             SL_DX(Vc,VX2,k,j - 1,i)/dx1(i)));
            #if DIMENSIONS == 3
              dVzi = SL::PLMLim(SL_DX(Vc,VX3,k,j,i + 1)/dx1(i+1),
                                  SL_DX(Vc,VX3,k,j - 1,i + 1)/dx1(i+1));
              dVzi = SL::PLMLim(dVzi,
                                  SL::PLMLim(SL_DX(Vc,VX3,k,j,i)/dx1(i),
                                               SL_DX(Vc,VX3,k,j - 1,i)/dx1(i)));
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
              dVxk = SL::PLMLim(SL_DZ(Vc,VX1,k + 1,j,i)/dx3(k+1),
                                  SL_DZ(Vc,VX1,k + 1,j - 1,i)/dx3(k+1));
              dVxk = SL::PLMLim(dVxk,
                                  SL::PLMLim(SL_DZ(Vc,VX1,k,j,i)/dx3(k),
                                               SL_DZ(Vc,VX1,k,j - 1,i)/dx3(k)));
              dVyk = SL::PLMLim(SL_DZ(Vc,VX2,k + 1,j,i)/dx3(k+1),
                                  SL_DZ(Vc,VX2,k + 1,j - 1,i)/dx3(k+1));
              dVyk = SL::PLMLim(dVyk,
                                  SL::PLMLim(SL_DZ(Vc,VX2,k,j,i)/dx3(k),
                                               SL_DZ(Vc,VX2,k,j - 1,i)/dx3(k)));
              dVzk = SL::PLMLim(SL_DZ(Vc,VX3,k + 1,j,i)/dx3(k+1),
                                  SL_DZ(Vc,VX3,k + 1,j - 1,i)/dx3(k+1));
              dVzk = SL::PLMLim(dVzk,
                                  SL::PLMLim(SL_DZ(Vc,VX3,k,j,i)/dx3(k),
                                               SL_DZ(Vc,VX3,k,j - 1,i)/dx3(k)));
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

          divV = D_EXPAND(dVxi, + dVyj, + dVzk);

          Pnor_par = 3.0*etaBrag*(bbgradV - divV/3.0);

          tau_xy = Pnor_par*bi*bj;
          tau_yy = Pnor_par*(bj*bj - 1./3.);
          tau_yz = Pnor_par*bj*bk;
        // No source term in cartesian geometry
        #else
          [[maybe_unused]] real vx1i, vx2i, vx3i;
          vx1i = vx2i = vx3i = ZERO_F;

          EXPAND(  vx1i = 0.5*(Vc(VX1,k,j-1,i)+Vc(VX1,k,j,i)); ,
                   vx2i = 0.5*(Vc(VX2,k,j-1,i)+Vc(VX2,k,j,i)); ,
                   vx3i = 0.5*(Vc(VX3,k,j-1,i)+Vc(VX3,k,j,i)); )

          Pnor_parC = 3.0*etaBragC*(bbgradVC - divVC/3.0);
        #endif // GEOMETRY != CARTESIAN

        #if GEOMETRY == CYLINDRICAL
          bbgradV = D_EXPAND(
                       bi*bi*dVxi + bi*bj*dVyi + bi*bk*(dVzi  - vx3i/x1(i)),
                     + bj*bi*dVxj + bj*bj*dVyj + bj*bk*dVzj + bk*bk*vx1i/x1(i),
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
                                                                                  + bk*bi*dVxk,
                    + bi*bj*(dVyi - vx2i/x1(i)) + bj*bj*(1./x1(i)*dVyj + vx1i/x1(i))
                                                                                  + bk*bj*dVyk,
                    + bi*bk*dVzi                         + bj*bk*1./x1(i)*dVzj
                                                                                  + bk*bk*dVzk
                    );

          divV = D_EXPAND(dVxi + vx1i/x1(i) ,
                           + 1./x1(i)*dVyj ,
                           + dVzk);

          // See Geoffroy's trick for curvature source terms in curvilinear coordinates
          EXPAND( bragViscSrc(IDIR,k,j,i) = ZERO_F;  ,
                  bragViscSrc(JDIR,k,j,i) = ZERO_F;  ,
                  bragViscSrc(KDIR,k,j,i) = ZERO_F;  )

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

          divV = D_EXPAND( 2.0*vx1i/x1(i) + dVxi,
                          +(SIN(x2(j))*Vc(VX2,k,j,i) - FABS(SIN(x2(j-1)))*Vc(VX2,k,j-1,i))/x1(i)
                           *one_dmu(j) ,
                          + dVzk/x1(i)*s_1 );

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
            etaBrag = 2.*(etaBragArr(k-1,j,i)*etaBragArr(k,j,i)) /
                         (etaBragArr(k-1,j,i)+etaBragArr(k,j,i));
          } else {
          etaBrag = HALF_F*(etaBragArr(k - 1,j,i)+etaBragArr(k,j,i));
          }
        } else {
          etaBragC = etaBrag = etaBragConstant;
        }

        if (haveSlopeLimiter) {
          dVxi = SL::PLMLim(SL_DX(Vc,VX1,k,j,i + 1)/dx1(i+1),
                              SL_DX(Vc,VX1,k - 1,j,i + 1)/dx1(i+1));
          dVxi = SL::PLMLim(dVxi,
                              SL::PLMLim(SL_DX(Vc,VX1,k,j,i)/dx1(i),
                                           SL_DX(Vc,VX1,k - 1,j,i)/dx1(i)));
          dVyi = SL::PLMLim(SL_DX(Vc,VX2,k,j,i + 1)/dx1(i+1),
                              SL_DX(Vc,VX2,k - 1,j,i + 1)/dx1(i+1));
          dVyi = SL::PLMLim(dVyi,
                              SL::PLMLim(SL_DX(Vc,VX2,k,j,i)/dx1(i),
                                           SL_DX(Vc,VX2,k - 1,j,i)/dx1(i)));
          dVzi = SL::PLMLim(SL_DX(Vc,VX3,k,j,i + 1)/dx1(i+1),
                              SL_DX(Vc,VX3,k - 1,j,i + 1)/dx1(i+1));
          dVzi = SL::PLMLim(dVzi,
                              SL::PLMLim(SL_DX(Vc,VX3,k,j,i)/dx1(i),
                                           SL_DX(Vc,VX3,k - 1,j,i)/dx1(i)));
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
            dVxj = SL::PLMLim(SL_DY(Vc,VX1,k,j + 1,i)/dx2(j+1),
                                SL_DY(Vc,VX1,k - 1,j + 1,i)/dx2(j+1));
            dVxj = SL::PLMLim(dVxj,
                                SL::PLMLim(SL_DY(Vc,VX1,k,j,i)/dx2(j),
                                             SL_DY(Vc,VX1,k - 1,j,i)/dx2(j)));
            dVyj = SL::PLMLim(SL_DY(Vc,VX2,k,j + 1,i)/dx2(j+1),
                                SL_DY(Vc,VX2,k - 1,j + 1,i)/dx2(j+1));
            dVyj = SL::PLMLim(dVyj,
                                SL::PLMLim(SL_DY(Vc,VX2,k,j,i)/dx2(j),
                                             SL_DY(Vc,VX2,k - 1,j,i)/dx2(j)));
            dVzj = SL::PLMLim(SL_DY(Vc,VX3,k,j + 1,i)/dx2(j+1),
                                SL_DY(Vc,VX3,k - 1,j + 1,i)/dx2(j+1));
            dVzj = SL::PLMLim(dVzj,
                                SL::PLMLim(SL_DY(Vc,VX3,k,j,i)/dx2(j),
                                             SL_DY(Vc,VX3,k - 1,j,i)/dx2(j)));
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

          divV = D_EXPAND(dVxi, + dVyj, + dVzk);

          Pnor_par = 3.0*etaBrag*(bbgradV - divV/3.0);

          tau_xz = Pnor_par*bi*bk;
          tau_yz = Pnor_par*bk*bj;
          tau_zz = Pnor_par*(bk*bk - 1./3.);

        // No source term in cartesian geometry
        #else
          [[maybe_unused]] real vx1i, vx2i, vx3i;
          vx1i = vx2i = vx3i = ZERO_F;

          EXPAND(  vx1i = 0.5*(Vc(VX1,k-1,j,i)+Vc(VX1,k,j,i)); ,
                   vx2i = 0.5*(Vc(VX2,k-1,j,i)+Vc(VX2,k,j,i)); ,
                   vx3i = 0.5*(Vc(VX3,k-1,j,i)+Vc(VX3,k,j,i)); )
        #endif // GEOMETRY != CARTESIAN

        #if GEOMETRY == CYLINDRICAL
          bbgradV = D_EXPAND(
                       bi*bi*dVxi + bi*bj*dVyi + bi*bk*(dVzi  - vx3i/x1(i)),
                     + bj*bi*dVxj + bj*bj*dVyj + bj*bk*dVzj + bk*bk*vx1i/x1(i),
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
                                                                                  + bk*bi*dVxk,
                    + bi*bj*(dVyi - vx2i/x1(i)) + bj*bj*(1./x1(i)*dVyj + vx1i/x1(i))
                                                                                  + bk*bj*dVyk,
                    + bi*bk*dVzi                         + bj*bk*1./x1(i)*dVzj
                                                                                  + bk*bk*dVzk
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


#undef D_DX_I
#undef D_DY_J
#undef D_DZ_K
#undef SL_DX
#undef SL_DY
#undef SL_DZ
#undef D_DY_I
#undef D_DZ_I
#undef D_DX_J
#undef D_DZ_J
#undef D_DX_K
#undef D_DY_K
#undef BX_I
#undef BX_J
#undef BX_K
#undef BY_I
#undef BY_J
#undef BY_K
#undef BZ_I
#undef BZ_J
#undef BZ_K
#endif // FLUID_BRAGINSKII_BRAGVISCOSITY_HPP_
