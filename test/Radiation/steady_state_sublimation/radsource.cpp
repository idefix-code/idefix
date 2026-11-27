// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file test/Radiation/steady_state_sublimation/radsource.cpp
//
// Last modified : 09/2026
//
// Copyright(C) by :
// - Nicolas Scepi <nicolas.scepi@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
// - Sebastien Valat <sebastien.valat@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "../idefix.hpp"
#include "radsource.hpp"
#include "physics.hpp"
#include "units.hpp"
#include "lookupTable.hpp"
#include "column.hpp"

void RadSource::RelativistCorrection(const real dt) {
  idfx::pushRegion("RadSource::RelativistCorrection");

  auto UcGas = this->UcGas;
  auto VcGas = this->VcGas;
  auto UcRad = this->UcRad;
  auto VcRad = this->VcRad;

  auto units = idfx::units;

  auto kp1D = this->kappa_planck_1D;
  auto kr1D = this->kappa_ross_1D;
  auto xi1D = this->xi_1D;

  auto kp2D = this->kappa_planck_2D;
  auto kr2D = this->kappa_ross_2D;
  auto xi2D = this->xi_2D;

  EquationOfState eos = this->eos;

  // Local copy of opacity parameters
  const Opacity kappa_type = this->kappa_type;
  const int kappa_ndim = this->kappa_ndim;
  const int xi_ndim = this->xi_ndim;

  IdefixArray3D<real> kappapArr;
  IdefixArray3D<real> kapparArr;
  IdefixArray3D<real> xiArr;
  real kappap_0,kappar_0,rho_0,T_0,xi_0;
  if (kappa_type == Opacity::constant) {
    kappap_0 = this->kappap_0;
    kappar_0 = this->kappar_0;
  } else if (kappa_type == Opacity::kramers) {
    kappap_0 = this->kappap_0;
    kappar_0 = this->kappar_0;
    T_0 = this->T_0;
    rho_0 = this->rho_0;
  } else if (kappa_type == Opacity::userfunc) {
    kappapArr = this->kappapArr;
    kapparArr = this->kapparArr;
  }

  const Opacity xi_type = this->xi_type;
  if (xi_type == Opacity::constant) {
    xi_0 = this->xi_0;
  } else if (xi_type == Opacity::userfunc) {
    xiArr = this->xiArr;
  }

  real reduced_c = this->reduced_c;


   idefix_for("RadSourceRelativistCorrection",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {

      real URad[RadiationPhysics::nvar];
      real UGas[DefaultPhysics::nvar];
      real VRad[RadiationPhysics::nvar];
      real VGas[DefaultPhysics::nvar];

      real kappa_p, kappa_r, xi;

      for(int nv = 0 ; nv < RadiationPhysics::nvar ; nv++) {
        URad[nv] = UcRad(nv,k,j,i);
        VRad[nv] = VcRad(nv,k,j,i);
      }
      for(int nv = 0 ; nv < DefaultPhysics::nvar ; nv++) {
        VGas[nv] = VcGas(nv,k,j,i);
        UGas[nv] = UcGas(nv,k,j,i);
      }

      // Compute total modified energy and momentum
      real Etot = UGas[ENG]+URad[ER]*units.c/(reduced_c*units.GetVelocity());
      EXPAND(real m1tot = UGas[MX1]+URad[FR1]/reduced_c;,
             real m2tot = UGas[MX2]+URad[FR2]/reduced_c;,
             real m3tot = UGas[MX3]+URad[FR3]/reduced_c;)

      real mu = eos.GetMu(VGas[PRS],VGas[RHO]);
      real T = VGas[PRS]/(VGas[RHO])*units.GetKelvin()*mu;

      // Compute opacities
      if (kappa_type == Opacity::constant) {
        kappa_p = kappap_0;
        kappa_r = kappar_0;
      } else if (kappa_type == Opacity::kramers) {
        kappa_p = kappap_0*VGas[RHO]*units.GetDensity()/rho_0*std::pow(T/T_0,-3.5);
        kappa_r = kappar_0*VGas[RHO]*units.GetDensity()/rho_0*std::pow(T/T_0,-3.5);
      } else if (kappa_type == Opacity::usertable) {
        real logT = std::log10(T);
        real logrho = std::log10(VGas[RHO]*units.GetDensity());
        if (kappa_ndim == 1) {
          kappa_p = kp1D.Get(&logT);
          kappa_r = kr1D.Get(&logT);
        } else if (kappa_ndim == 2) {
          real x[2];
          x[1] = FMIN(FMAX(logT,2.5),5.98);
          x[0] = FMIN(-4.05,FMAX(-14.,logrho));
          kappa_p = std::pow(10.,kp2D.Get(x));
          kappa_r = std::pow(10.,kr2D.Get(x));
        }
      } else if (kappa_type == Opacity::userfunc) {
        kappa_p = kappapArr(k,j,i);
        kappa_r = kapparArr(k,j,i);
      }

      if (xi_type == Opacity::constant) {
        xi = xi_0;
      } else if (xi_type == Opacity::usertable) {
        real logT = std::log10(T);
        real logrho = std::log10(VGas[RHO]*units.GetDensity());
        if (xi_ndim == 1) {
          xi = xi1D.Get(&logT);
        } else if (xi_ndim == 2) {
          real x[2];
          x[1] = logT;
          x[0] = logrho;
          xi = std::pow(10.,xi2D.Get(x));
        }
      } else if (xi_type== Opacity::userfunc) {
        xi = xiArr(k,j,i);
      }

      // Compute beta parameter vector
      EXPAND(real beta1 = VGas[VX1]*units.GetVelocity()/units.c;,
             real beta2 = VGas[VX2]*units.GetVelocity()/units.c;,
             real beta3 = VGas[VX3]*units.GetVelocity()/units.c;)

      // Dot product of beta and Fr
      real betaFr = EXPAND(VRad[FR1]*beta1,+VRad[FR2]*beta2,+VRad[FR3]*beta3);

      // Compute beta**2
      real betasq = EXPAND(beta1*beta1,+beta2*beta2,+beta3*beta3);

      // Compute radiation pressure tensor
      real Fnorm2 = EXPAND(VRad[FR1]*VRad[FR1] , + VRad[FR2]*VRad[FR2], + VRad[FR3]*VRad[FR3]);
      real inv_Fnorm2 = (Fnorm2 <= 1.e-100 ? 1.e-100 : ONE_F / Fnorm2);
      real Er2 = VRad[ER]*VRad[ER];
      real f_param2 = (Er2 < 1.e-100 ? Fnorm2/(1.e-100) : Fnorm2/(Er2));
      real chi  = 3.+4.*f_param2;
      chi /= 5.+2.*std::sqrt(4.-3.*f_param2);

      // Add momentum-like part of the radiation pressure tensor
      EXPAND ( real P11 = HALF_F*(3.*chi-1.)*VRad[ER]*VRad[FR1]*VRad[FR1]*inv_Fnorm2;,
               real P12 = HALF_F*(3.*chi-1.)*VRad[ER]*VRad[FR1]*VRad[FR2]*inv_Fnorm2;
               real P22 = HALF_F*(3.*chi-1.)*VRad[ER]*VRad[FR2]*VRad[FR2]*inv_Fnorm2;,
               real P13 = HALF_F*(3.*chi-1.)*VRad[ER]*VRad[FR1]*VRad[FR3]*inv_Fnorm2;
               real P23 = HALF_F*(3.*chi-1.)*VRad[ER]*VRad[FR2]*VRad[FR3]*inv_Fnorm2;
               real P33 = HALF_F*(3.*chi-1.)*VRad[ER]*VRad[FR3]*VRad[FR3]*inv_Fnorm2;)

      // Add pressure-like part of the radiation pressure tensor
      EXPAND ( P11 += HALF_F*(1.-chi)*VRad[ER];,
               P22 += HALF_F*(1.-chi)*VRad[ER];,
               P33 += HALF_F*(1.-chi)*VRad[ER];)

      // Compute beta.(beta.P)
      real beta2P = EXPAND(beta1*beta1*P11, +2.*beta1*beta2*P12+beta2*beta2*P22, +2.*beta1*beta3*P13+2.*beta2*beta3*P23+beta3*beta3*P33);

      // Add relativistic correction to energy source term
      real G0 = -2.*betaFr*kappa_p;
      G0 += (xi+kappa_r)*(betaFr - betasq*VRad[ER] - beta2P);
      G0 *= VGas[RHO]*units.GetDensity();

      // Add relativistic correction to flux source term
      EXPAND ( real G1 = kappa_p*beta1*(VRad[ER]-units.ar*std::pow(T,4)/units.GetEnergy()-2.*betaFr);,
              real G2 = kappa_p*beta2*(VRad[ER]-units.ar*std::pow(T,4)/units.GetEnergy()-2.*betaFr);,
               real G3 = kappa_p*beta3*(VRad[ER]-units.ar*std::pow(T,4)/units.GetEnergy()-2.*betaFr);)


      EXPAND ( G1 -= (xi+kappa_r)*(beta1*P11+VRad[ER]*beta1);,
               G1 -= (xi+kappa_r)*beta2*P12;
               G2 -= (xi+kappa_r)*(beta1*P12+beta2*P22+VRad[ER]*beta2);,
               G1 -= (xi+kappa_r)*beta3*P13;
               G2 -= (xi+kappa_r)*beta3*P23;
               G3 -= (xi+kappa_r)*(beta1*P13+beta2*P23+beta3*P33+VRad[ER]*beta3);)
      EXPAND ( G1 *= VGas[RHO]*units.GetDensity();,
               G2 *= VGas[RHO]*units.GetDensity();,
               G3 *= VGas[RHO]*units.GetDensity();)


      URad[ER] -= G0*dt*units.GetTime()*reduced_c*units.GetVelocity();
      EXPAND( URad[FR1] -= G1*dt*units.GetTime()*reduced_c*units.GetVelocity();,
              URad[FR2] -= G2*dt*units.GetTime()*reduced_c*units.GetVelocity();,
              URad[FR3] -= G3*dt*units.GetTime()*reduced_c*units.GetVelocity();)

      if ((Etot - URad[ER]*units.c/(reduced_c*units.GetVelocity()))<=ZERO_F) {
        #ifdef SMALL_ER
        URad[ER] = SMALL_ER;
        //std::printf("URad[ER]=%e UGas[ENG]=%e Etot=%e T=%e rho=%e kappap=%e at i=%i j=%i and k=%i\n",URad[ER],UGas[ENG],Etot,T,VGas[RHO]*units.GetDensity(),kappa_p,i,j,k);
        //UGas[ENG] = Etot - URad[ER]*units.c/(reduced_c*units.GetVelocity());
        #else
        Kokkos::abort("ENG=0 in RadSourceRelativistCorrection");
        #endif
      } else {
        UGas[ENG] = Etot - URad[ER]*units.c/(reduced_c*units.GetVelocity());
      }
      EXPAND( UGas[MX1] = m1tot - URad[FR1]/reduced_c;,
              UGas[MX2] = m2tot - URad[FR2]/reduced_c;,
              UGas[MX3] = m3tot - URad[FR3]/reduced_c;)

      for(int nv = 0 ; nv < RadiationPhysics::nvar ; nv++) {
        UcRad(nv,k,j,i) = URad[nv];
      }

      for(int nv = 0 ; nv < DefaultPhysics::nvar ; nv++) {
        UcGas(nv,k,j,i) = UGas[nv];
      }

  });



}

void RadSource::SourceFullImplicit(const real dt) {
  idfx::pushRegion("RadSource::Source_full_implicit");

  // Ensure that radiation cannot be run with isothermal eos
  #ifndef ISOTHERMAL

  auto UcGas = this->UcGas;
  auto VcGas = this->VcGas;
  auto UcRad = this->UcRad;
  auto VcRad = this->VcRad;
  auto InvDt = this->InvDt;

  auto units = idfx::units;

  auto kp1D = this->kappa_planck_1D;
  auto kr1D = this->kappa_ross_1D;
  auto xi1D = this->xi_1D;

  auto kp2D = this->kappa_planck_2D;
  auto kr2D = this->kappa_ross_2D;
  auto xi2D = this->xi_2D;

  EquationOfState eos = this->eos;

  real reduced_c = this->reduced_c;

  // Irradiation source
  bool irr_flag=false;
  IdefixArray3D<real> divF = this->divF;
  if (haveIrradiation){
    IrrFlux(divF);
    irr_flag=true;
  }

  IdefixArray3D<real> Qvisc = this->Qvisc;
  Qviscous(Qvisc);


  // Local copy of opacity parameters
  const Opacity kappa_type = this->kappa_type;
  const int kappa_ndim = this->kappa_ndim;
  const int xi_ndim = this->xi_ndim;

  IdefixArray3D<real> kappapArr;
  IdefixArray3D<real> kapparArr;
  IdefixArray3D<real> xiArr;
  real kappap_0,kappar_0,rho_0,T_0,xi_0;
  if (kappa_type == Opacity::constant) {
    kappap_0 = this->kappap_0;
    kappar_0 = this->kappar_0;
  } else if (kappa_type == Opacity::kramers) {
    kappap_0 = this->kappap_0;
    kappar_0 = this->kappar_0;
    T_0 = this->T_0;
    rho_0 = this->rho_0;
  } else if (kappa_type == Opacity::userfunc) {
    kappapArr = this->kappapArr;
    kapparArr = this->kapparArr;
  }

  const Opacity xi_type = this->xi_type;
  if (xi_type == Opacity::constant) {
    xi_0 = this->xi_0;
  } else if (xi_type == Opacity::userfunc) {
    xiArr = this->xiArr;
  }


  idefix_for("RadSourceFullImplicit",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {

      real URad[RadiationPhysics::nvar];
      real UGas[DefaultPhysics::nvar];
      real VGas[DefaultPhysics::nvar];

      real kappa_p, kappa_r, xi;

      for(int nv = 0 ; nv < RadiationPhysics::nvar ; nv++) {
        URad[nv] = UcRad(nv,k,j,i);
      }
      for(int nv = 0 ; nv < DefaultPhysics::nvar ; nv++) {
        UGas[nv] = UcGas(nv,k,j,i);
        VGas[nv] = VcGas(nv,k,j,i);
      }

      // Compute total modified energy
      real Etot = UGas[ENG]+URad[ER]*units.c/(reduced_c*units.GetVelocity());

      // Add irradiation heating if needed
      if (irr_flag){
        Etot -= divF(k,j,i)*dt*units.GetTime()/units.GetEnergy();
      }

      // Add viscous energy
      Etot += Qvisc(k,j,i)*dt;

      // Compute total modified momentum
      EXPAND(real m1tot = UGas[MX1]+URad[FR1]/reduced_c;,
             real m2tot = UGas[MX2]+URad[FR2]/reduced_c;,
             real m3tot = UGas[MX3]+URad[FR3]/reduced_c;)

      // Store conserved variables after hyperbolic step
      real Er_hyp = URad[ER];
      EXPAND(real Fr1_hyp = URad[FR1];,
             real Fr2_hyp = URad[FR2];,
             real Fr3_hyp = URad[FR3];)

      real Fnorm = std::sqrt(EXPAND(URad[FR1]*URad[FR1] , + URad[FR2]*URad[FR2], + URad[FR3]*URad[FR3]));

      real mu = eos.GetMu(VGas[PRS],VGas[RHO]);
      real cv = units.k_B/(units.u*mu);

      real T = VGas[PRS]/(VGas[RHO])*units.GetKelvin()*mu;
      real T3 = std::pow(T,3);

      // Compute opacities
      if (kappa_type == Opacity::constant) {
        kappa_p = kappap_0;
        kappa_r = kappar_0;
      } else if (kappa_type == Opacity::kramers) {
        kappa_p = kappap_0*VGas[RHO]*units.GetDensity()/rho_0*std::pow(T/T_0,-3.5);
        kappa_r = kappar_0*VGas[RHO]*units.GetDensity()/rho_0*std::pow(T/T_0,-3.5);
      } else if (kappa_type == Opacity::usertable) {
        real logT = std::log10(T);
        real logrho = std::log10(VGas[RHO]*units.GetDensity());
        if (kappa_ndim == 1) {
          kappa_p = kp1D.Get(&logT);
          kappa_r = kr1D.Get(&logT);
        } else if (kappa_ndim == 2) {
          real x[2];
          x[1] = FMIN(FMAX(logT,2.5),5.98);
          x[0] = FMIN(-4.05,FMAX(-14.,logrho));
          kappa_p = std::pow(10.,kp2D.Get(x));
          kappa_r = std::pow(10.,kr2D.Get(x));
        }
      } else if (kappa_type == Opacity::userfunc) {
        kappa_p = kappapArr(k,j,i);
        kappa_r = kapparArr(k,j,i);
      }

      if (xi_type == Opacity::constant) {
        xi = xi_0;
      } else if (xi_type == Opacity::usertable) {
        real logT = std::log10(T);
        real logrho = std::log10(VGas[RHO]*units.GetDensity());
        if (xi_ndim == 1) {
          xi = xi1D.Get(&logT);
        } else if (xi_ndim == 2) {
          real x[2];
          x[1] = logT;
          x[0] = logrho;
          xi = std::pow(10.,xi2D.Get(x));
        }
      } else if (xi_type== Opacity::userfunc) {
        xi = xiArr(k,j,i);
      }


      real kk_red = reduced_c * units.GetVelocity() * dt * units.GetTime() * kappa_p * VGas[RHO]*units.GetDensity();
      real kk = units.c * dt * units.GetTime() * kappa_p * VGas[RHO]*units.GetDensity();
      real xx_red = reduced_c * units.GetVelocity() * dt * units.GetTime() * (xi + kappa_r) * VGas[RHO]*units.GetDensity();

      // Define matrix to invert
      real gamma = eos.GetGamma(VGas[PRS],VGas[RHO]);

      real M00 = ONE_F + kk_red;
      real M11 = VGas[RHO]*units.GetDensity()*cv/(gamma-1.) + 4.*kk*units.ar*T3;
      real M01 = -4.*kk_red*units.ar*T3;
      real M10 = -kk;

      // Define right-hand side of system
      real S0 = Er_hyp*units.GetEnergy() - 3.*kk_red*units.ar*T3*T;
      real S1 = VGas[RHO]*units.GetDensity()*cv*T/(gamma-1.) + 3.*kk*units.ar*T3*T;

      // Add irradiation heating to RHS if needed
      if (irr_flag){
        S1 -= divF(k,j,i)*dt*units.GetTime();
      }

      // Add viscous heating to RHS
      S1 += Qvisc(k,j,i)*dt*units.GetEnergy();

      // Invert system
      real det = M00*M11 - M01*M10;

      real Minv00 = M11/det;
      //real Minv11 = M00/det;
      real Minv01 = -M01/det;
      //real Minv10 = -M10/det;

      real Er_new = Minv00*S0 + Minv01*S1;
      //real T_new = Minv10*S0 + Minv11*S1;

      // Update conservative variables
      URad[ER] = Er_new/units.GetEnergy();
      EXPAND( URad[FR1] = Fr1_hyp/(1.+xx_red);,
              URad[FR2] = Fr2_hyp/(1.+xx_red);,
              URad[FR3] = Fr3_hyp/(1.+xx_red);)
      //if (i==2 && j==495) std::printf("Etot -ERad*c/redc=%e URad[ER]=%e Er_hyp=%e UGas[ENG]=%e Etot=%e T=%e arT4=%e rho=%e kappap=%e kappar=%e at i=%i j=%i and k=%i\n",(Etot - URad[ER]*units.c/(reduced_c*units.GetVelocity())),URad[ER],Er_hyp,UGas[ENG],Etot,T,units.ar*T3*T/units.GetEnergy(),VGas[RHO]*units.GetDensity(),kappa_p,kappa_r,i,j,k);

      if ((Etot - URad[ER]*units.c/(reduced_c*units.GetVelocity()))<=ZERO_F) {
        #ifdef SMALL_ER
        //std::printf("URad[ER]=%e Er_hyp=%e UGas[ENG]=%e Etot=%e T=%e arT4=%e rho=%e kappap=%e kappar=%e at i=%i j=%i and k=%i\n",URad[ER],Er_hyp,UGas[ENG],Etot,T,units.ar*T3*T/units.GetEnergy(),VGas[RHO]*units.GetDensity(),kappa_p,kappa_r,i,j,k);
        URad[ER] = units.ar*T3*T/units.GetEnergy();
        UGas[ENG] = Etot - URad[ER]*units.c/(reduced_c*units.GetVelocity());
        #else
        Kokkos::abort("ENG=0 in RadSourceFullImplicit");
        #endif

      } else {
        UGas[ENG] = Etot - URad[ER]*units.c/(reduced_c*units.GetVelocity());
      }

      EXPAND( UGas[MX1] = m1tot - URad[FR1]/reduced_c;,
              UGas[MX2] = m2tot - URad[FR2]/reduced_c;,
              UGas[MX3] = m3tot - URad[FR3]/reduced_c;)

      for(int nv = 0 ; nv < RadiationPhysics::nvar ; nv++) {
        UcRad(nv,k,j,i) = URad[nv];
      }

      for(int nv = 0 ; nv < DefaultPhysics::nvar ; nv++) {
        UcGas(nv,k,j,i) = UGas[nv];
      }
    });

    #endif
  idfx::popRegion();
}


void RadSource::SourceFixedPointRad(const real dt) {
  idfx::pushRegion("RadSource::SourceFixedPointRad");

  // Ensure that radiation cannot be run with isothermal eos
  #ifndef ISOTHERMAL

  auto UcGas = this->UcGas;
  auto VcGas = this->VcGas;
  auto UcRad = this->UcRad;
  auto VcRad = this->VcRad;
  auto InvDt = this->InvDt;

  auto units=idfx::units;

  auto kp1D = this->kappa_planck_1D;
  auto kr1D = this->kappa_ross_1D;
  auto xi1D = this->xi_1D;

  auto kp2D = this->kappa_planck_2D;
  auto kr2D = this->kappa_ross_2D;
  auto xi2D = this->xi_2D;

  real reduced_c = this->reduced_c;

  // Max iteration for fixed-point solver
  int MAX_ITER = 200;
  // Tolerance on ER and ENG for fixed-point solver
  real tol = 1.e-3;

  EquationOfState eos = this->eos;

  // Irradiation source
  bool irr_flag=false;
  IdefixArray3D<real> divF = this->divF;
  if (haveIrradiation){
    IrrFlux(divF);
    irr_flag=true;
  }

  // Local copy of opacity parameters
  const Opacity kappa_type = this->kappa_type;
  const int kappa_ndim = this->kappa_ndim;
  const int xi_ndim = this->xi_ndim;

  real kappap_0,kappar_0,rho_0,T_0,xi_0;
  IdefixArray3D<real> kappapArr;
  IdefixArray3D<real> kapparArr;
  IdefixArray3D<real> xiArr;
  if (kappa_type == Opacity::constant) {
    kappap_0 = this->kappap_0;
    kappar_0 = this->kappar_0;
  } else if (kappa_type == Opacity::kramers) {
    kappap_0 = this->kappap_0;
    kappar_0 = this->kappar_0;
    T_0 = this->T_0;
    rho_0 = this->rho_0;
  } else if (kappa_type == Opacity::userfunc) {
    kappapArr = this->kappapArr;
    kapparArr = this->kapparArr;
  }

  const Opacity xi_type = this->xi_type;
  if (xi_type == Opacity::constant) {
    xi_0 = this->xi_0;
  } else if (xi_type == Opacity::userfunc) {
    xiArr = this->xiArr;
  }

  idefix_for("RadSourceFixedPointRad",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {

      // Add heating due to irradiation flux if needed
      if (irr_flag){
        // Limit time step relative to characteristic time of irradiation heating
        InvDt(k,j,i) += FABS(units.GetTime()*divF(k,j,i)/units.GetEnergy()/UcGas(ENG,k,j,i));
        UcGas(ENG,k,j,i) -= dt*units.GetTime()*divF(k,j,i)/units.GetEnergy();
      }

      real UGas[DefaultPhysics::nvar];
      real VGas[DefaultPhysics::nvar];
      real URad[RadiationPhysics::nvar];

      real kappa_p, kappa_r, xi;

      for(int nv = 0 ; nv < RadiationPhysics::nvar ; nv++) {
        URad[nv] = UcRad(nv,k,j,i);
      }

      for(int nv = 0 ; nv < DefaultPhysics::nvar ; nv++) {
        UGas[nv] = UcGas(nv,k,j,i);
        VGas[nv] = VcGas(nv,k,j,i);
      }

      // Compute total modified energy and momentum
      real Etot = UGas[ENG]+URad[ER]*units.c/(reduced_c*units.GetVelocity());
      EXPAND(real m1tot = UGas[MX1]+URad[FR1]/reduced_c;,
             real m2tot = UGas[MX2]+URad[FR2]/reduced_c;,
             real m3tot = UGas[MX3]+URad[FR3]/reduced_c;)

      // Store conservative variables after hydro step
      real Er_hyp = URad[ER];
      EXPAND(real Fr1_hyp = URad[FR1];,
             real Fr2_hyp = URad[FR2];,
             real Fr3_hyp = URad[FR3];)

      real Er_old, Fnorm_old, Egas_old, Mnorm_old;

      real err1= 1.;
      real err2= 1.;
      real err3= 1.;
      real err4= 1.;
      int count = 0;

      real Fnorm = std::sqrt(EXPAND(URad[FR1]*URad[FR1] , + URad[FR2]*URad[FR2], + URad[FR3]*URad[FR3]));
      real Mnorm = std::sqrt(EXPAND(UGas[MX1]*UGas[MX1] , + UGas[MX2]*UGas[MX2], + UGas[MX3]*UGas[MX3]));

      // Assume mu is constant during iteration (to check)
      real mu = eos.GetMu(VGas[PRS],VGas[RHO]);

      real T = VGas[PRS]/(VGas[RHO])*units.GetKelvin()*mu;

      // Assume kappa and xi are constant during iteration (to check)
      if (kappa_type == Opacity::constant) {
        kappa_p = kappap_0;
        kappa_r = kappar_0;
      } else if (kappa_type == Opacity::kramers) {
        kappa_p = kappap_0*VGas[RHO]*units.GetDensity()/rho_0*std::pow(T/T_0,-3.5);
        kappa_r = kappar_0*VGas[RHO]*units.GetDensity()/rho_0*std::pow(T/T_0,-3.5);
      } else if (kappa_type == Opacity::usertable) {
        real logT = std::log10(T);
        real logrho = std::log10(VGas[RHO]*units.GetDensity());
        if (kappa_ndim == 1) {
          kappa_p = kp1D.Get(&logT);
          kappa_r = kr1D.Get(&logT);
        } else if (kappa_ndim == 2) {
          real x[2];
          x[1] = FMIN(FMAX(logT,2.5),5.98);
          x[0] = FMIN(-4.05,FMAX(-14.,logrho));
          kappa_p = std::pow(10.,kp2D.Get(x));
          kappa_r = std::pow(10.,kr2D.Get(x));
        }
      } else if (kappa_type == Opacity::userfunc) {
        kappa_p = kappapArr(k,j,i);
        kappa_r = kapparArr(k,j,i);
      }

      if (xi_type == Opacity::constant) {
        xi = xi_0;
      } else if (xi_type == Opacity::usertable) {
        real logT = std::log10(T);
        real logrho = std::log10(VGas[RHO]*units.GetDensity());
        if (xi_ndim == 1) {
          xi = xi1D.Get(&logT);
        } else if (xi_ndim == 2) {
          real x[2];
          x[1] = logT;
          x[0] = logrho;
          xi = std::pow(10.,xi2D.Get(x));
        }
      } else if (xi_type== Opacity::userfunc) {
        xi = xiArr(k,j,i);
      }

      // Iterate on radiative variables
      while (((err1>tol) || (err2>tol) || (err3>tol) || (err4>tol)) && (count < MAX_ITER)){

        Er_old = URad[ER];
        Fnorm_old = Fnorm;
        Egas_old = UGas[ENG];
        Mnorm_old = Mnorm;

        real kk_red = reduced_c * units.GetVelocity() * dt * units.GetTime() * kappa_p * VGas[RHO]*units.GetDensity();
        real xx_red = reduced_c * units.GetVelocity() * dt * units.GetTime() * (xi + kappa_r) * VGas[RHO]*units.GetDensity();

        // "Implicit" step on radiation conservative variables
        URad[ER] = Er_hyp +  kk_red*units.ar*std::pow(T,4.)/units.GetEnergy();
        URad[ER] /= 1. + kk_red;
        EXPAND( URad[FR1] = Fr1_hyp/(1.+xx_red);,
                URad[FR2] = Fr2_hyp/(1.+xx_red);,
                URad[FR3] = Fr3_hyp/(1.+xx_red);)
        Fnorm = std::sqrt(EXPAND(URad[FR1]*URad[FR1] , + URad[FR2]*URad[FR2], + URad[FR3]*URad[FR3]));

        // Update gas conservative variables
        if ((Etot - URad[ER]*units.c/(reduced_c*units.GetVelocity()))<=ZERO_F) {
          Kokkos::abort("ENG=0 in RadSourceFixedPointRad");
        } else {
          UGas[ENG] = Etot - URad[ER]*units.c/(reduced_c*units.GetVelocity());
        }
        EXPAND( UGas[MX1] = m1tot - URad[FR1]/reduced_c;,
                UGas[MX2] = m2tot - URad[FR2]/reduced_c;,
                UGas[MX3] = m3tot - URad[FR3]/reduced_c;)
        Mnorm = std::sqrt(EXPAND(UGas[MX1]*UGas[MX1] , + UGas[MX2]*UGas[MX2], + UGas[MX3]*UGas[MX3]));

        // Update gas primitive variables
        K_ConsToPrim<DefaultPhysics>(VGas, UGas, &eos);

        // Compute new temperature
        T = VGas[PRS]/(VGas[RHO])*units.GetKelvin()*mu;

        // Compute errors and number of cycles
        err1 = std::abs(1.-URad[ER]/Er_old);
        err2 = std::abs(1.-Fnorm/Fnorm_old);
        err3 = std::abs(1.-UGas[ENG]/Egas_old);
        err4 = std::abs(1.-Mnorm/Mnorm_old);
        count += 1;
      }

      for(int nv = 0 ; nv < RadiationPhysics::nvar ; nv++) {
        UcRad(nv,k,j,i) = URad[nv];
      }

      for(int nv = 0 ; nv < DefaultPhysics::nvar ; nv++) {
        UcGas(nv,k,j,i) = UGas[nv];
      }
  });

  #endif

  idfx::popRegion();
}


void RadSource::SourceFixedPointGas(const real dt) {
  idfx::pushRegion("RadSource::SourceFixedPointGas");

  // Ensure that radiation cannot be run with isothermal eos
  #ifndef ISOTHERMAL

  auto UcGas = this->UcGas;
  auto VcGas = this->VcGas;
  auto UcRad = this->UcRad;
  auto VcRad = this->VcRad;
  auto InvDt = this->InvDt;

  real reduced_c = this->reduced_c;

  auto units=idfx::units;

  auto kp1D = this->kappa_planck_1D;
  auto kr1D = this->kappa_ross_1D;
  auto xi1D = this->xi_1D;

  auto kp2D = this->kappa_planck_2D;
  auto kr2D = this->kappa_ross_2D;
  auto xi2D = this->xi_2D;

  // Max iteration for fixed-point solver
  int MAX_ITER = 200;
  // Tolerance on ER and ENG for fixed-point solver
  real tol = 1.e-3;

  EquationOfState eos = this->eos;

  // Local copy of opacity parameters
  const Opacity kappa_type = this->kappa_type;
  const int kappa_ndim = this->kappa_ndim;
  const int xi_ndim = this->xi_ndim;

  real kappap_0,kappar_0,rho_0,T_0,xi_0;
  IdefixArray3D<real> kappapArr;
  IdefixArray3D<real> kapparArr;
  IdefixArray3D<real> xiArr;
  if (kappa_type == Opacity::constant) {
    kappap_0 = this->kappap_0;
    kappar_0 = this->kappar_0;
  } else if (kappa_type == Opacity::kramers) {
    kappap_0 = this->kappap_0;
    kappar_0 = this->kappar_0;
    T_0 = this->T_0;
    rho_0 = this->rho_0;
  } else if (kappa_type == Opacity::userfunc) {
    kappapArr = this->kappapArr;
    kapparArr = this->kapparArr;
  }

  const Opacity xi_type = this->xi_type;
  if (xi_type == Opacity::constant) {
    xi_0 = this->xi_0;
  } else if (xi_type == Opacity::userfunc) {
    xiArr = this->xiArr;
  }

  idefix_for("RadSource",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {

      real Etot = UcGas(ENG,k,j,i)+UcRad(ER,k,j,i)*units.c/(reduced_c*units.GetVelocity());
      real m1tot = UcGas(MX1,k,j,i)+UcRad(FR1,k,j,i)/reduced_c;
      real m2tot = UcGas(MX2,k,j,i)+UcRad(FR2,k,j,i)/reduced_c;
      real m3tot = UcGas(MX3,k,j,i)+UcRad(FR3,k,j,i)/reduced_c;

      real Er_hyp = UcRad(ER,k,j,i);
      real Fr1_hyp = UcRad(FR1,k,j,i);
      real Fr2_hyp = UcRad(FR2,k,j,i);
      real Fr3_hyp = UcRad(FR3,k,j,i);

      real Egas_hyp = UcGas(ENG,k,j,i);
      real m1gas_hyp = UcGas(MX1,k,j,i);
      real m2gas_hyp = UcGas(MX2,k,j,i);
      real m3gas_hyp = UcGas(MX3,k,j,i);

      real URad[RadiationPhysics::nvar];
      real Er_old, Fnorm_old, Egas_old, Mnorm_old;

      real UGas[DefaultPhysics::nvar];
      real VGas[DefaultPhysics::nvar];

      real kappa_p, kappa_r, xi;

      for(int nv = 0 ; nv < RadiationPhysics::nvar ; nv++) {
        URad[nv] = UcRad(nv,k,j,i);
      }

      for(int nv = 0 ; nv < DefaultPhysics::nvar ; nv++) {
        UGas[nv] = UcGas(nv,k,j,i);
        VGas[nv] = VcGas(nv,k,j,i);
      }

      real err1= 1.;
      real err2= 1.;
      real err3= 1.;
      real err4= 1.;
      int count = 0;

      real Fnorm = std::sqrt(EXPAND(URad[FR1]*URad[FR1] , + URad[FR2]*URad[FR2], + URad[FR3]*URad[FR3]));
      real Mnorm = std::sqrt(EXPAND(UGas[MX1]*UGas[MX1] , + UGas[MX2]*UGas[MX2], + UGas[MX3]*UGas[MX3]));

      // Assume mu is constant during iteration (to check)
      real mu = eos.GetMu(VGas[PRS],VGas[RHO]);
      real T = VGas[PRS]/(VGas[RHO])*units.GetKelvin()*mu;

      // Compute opacities (out of while loop so that opacity is contant throughout the fixed_point iteration)
      if (kappa_type == Opacity::constant) {
        kappa_p = kappap_0;
        kappa_r = kappar_0;
      } else if (kappa_type == Opacity::kramers) {
        kappa_p = kappap_0*VGas[RHO]*units.GetDensity()/rho_0*std::pow(T/T_0,-3.5);
        kappa_r = kappar_0*VGas[RHO]*units.GetDensity()/rho_0*std::pow(T/T_0,-3.5);
      } else if (kappa_type == Opacity::usertable) {
        real logT = std::log10(T);
        real logrho = std::log10(VGas[RHO]*units.GetDensity());
        if (kappa_ndim == 1) {
          kappa_p = kp1D.Get(&logT);
          kappa_r = kr1D.Get(&logT);
        } else if (kappa_ndim == 2) {
          real x[2];
          x[1] = logT;
          x[0] = logrho;
          kappa_p = std::pow(10.,kp2D.Get(x));
          kappa_r = std::pow(10.,kr2D.Get(x));
        }
      } else if (kappa_type == Opacity::userfunc) {
        kappa_p = kappapArr(k,j,i);
        kappa_r = kapparArr(k,j,i);
      }

      if (xi_type == Opacity::constant) {
        xi = xi_0;
      } else if (xi_type == Opacity::usertable) {
        real logT = std::log10(T);
        real logrho = std::log10(VGas[RHO]*units.GetDensity());
        if (xi_ndim == 1) {
          xi = xi1D.Get(&logT);
        } else if (xi_ndim == 2) {
          real x[2];
          x[1] = logT;
          x[0] = logrho;
          xi = std::pow(10.,xi2D.Get(x));
        }
      } else if (xi_type== Opacity::userfunc) {
        xi = xiArr(k,j,i);
      }

      while (((err1>tol) || (err2>tol) || (err3>tol) || (err4>tol)) && (count < MAX_ITER)){

        Er_old = URad[ER];
        Fnorm_old = Fnorm;
        Egas_old = UGas[ENG];
        Mnorm_old = Mnorm;

        real kk =  units.c * dt * units.GetTime() * kappa_p * VGas[RHO]*units.GetDensity();
        real xx =  dt * units.GetTime() * (xi + kappa_r) * VGas[RHO]*units.GetDensity();

        // Stop if UGas <= 0
        if ((Egas_hyp +  kk*(URad[ER]-units.ar*std::pow(T,4)/units.GetEnergy()))<=ZERO_F) {
          Kokkos::abort("EGas=0 in Radsource");
          UGas[ENG] = 1.e-6;
        } else {
          UGas[ENG] = Egas_hyp +  kk*(URad[ER]-units.ar*std::pow(T,4)/units.GetEnergy());
        }

        EXPAND( UGas[MX1] = m1gas_hyp + URad[FR1]*xx;,
                UGas[MX2] = m2gas_hyp + URad[FR2]*xx;,
                UGas[MX3] = m3gas_hyp + URad[FR3]*xx;)
        Mnorm = std::sqrt(EXPAND(UGas[MX1]*UGas[MX1] , + UGas[MX2]*UGas[MX2], + UGas[MX3]*UGas[MX3]));

        URad[ER] = (Etot - UGas[ENG])*reduced_c*units.GetVelocity()/units.c;

        // Stop if URad <= 0
        if (URad[ER]<=ZERO_F) {
          Kokkos::abort("ERad=0 in RadSourceFixedPointGas");
        }

        EXPAND( URad[FR1] = (m1tot - UGas[MX1])*reduced_c;,
                URad[FR2] = (m2tot - UGas[MX2])*reduced_c;,
                URad[FR3] = (m3tot - UGas[MX3])*reduced_c;)

        Fnorm = std::sqrt(EXPAND(URad[FR1]*URad[FR1] , + URad[FR2]*URad[FR2], + URad[FR3]*URad[FR3]));

        K_ConsToPrim<DefaultPhysics>(VGas, UGas, &eos);
        T = VGas[PRS]/(VGas[RHO])*units.GetKelvin()*mu;

        err1 = std::abs(1.-UGas[ENG]/Egas_old);
        err2 = std::abs(1.-Mnorm/Mnorm_old);
        err3 = std::abs(1.-URad[ER]/Er_old);
        err4 = std::abs(1.-Fnorm/Fnorm_old);
        count += 1;

      }

      for(int nv = 0 ; nv < RadiationPhysics::nvar ; nv++) {
        UcRad(nv,k,j,i) = URad[nv];
      }

      for(int nv = 0 ; nv < DefaultPhysics::nvar ; nv++) {
        UcGas(nv,k,j,i) = UGas[nv];
      }
  });

  #endif

  idfx::popRegion();
}

void RadSource::ShowConfig() {

  // Ensure that radiation cannot be run with isothermal eos
  #ifndef ISOTHERMAL

  idfx::cout << "RadSource: kappa is ";
  switch(kappa_type) {
    case Opacity::constant:
      idfx::cout << "constant." << std::endl;
      break;
    case Opacity::kramers:
      idfx::cout << "from kramers' law." << std::endl;
      break;
    case Opacity::usertable:
      idfx::cout << "from a user table." << std::endl;
      break;
    case Opacity::userfunc:
      idfx::cout << "from a user-defined function."
                     << std::endl;
      if(!data->radiation[0]->kappaFunc) {
        IDEFIX_ERROR("No opacity function has been enrolled for kappa");
      }
      break;
  }
  idfx::cout << "RadSource: xi is ";
  switch(xi_type) {
    case Opacity::constant:
      idfx::cout << "constant." << std::endl;
      break;
    case Opacity::kramers:
      idfx::cout << "!!! from kramers' law, which is not allowed !!!" << std::endl;
      break;
    case Opacity::usertable:
      idfx::cout << "from a user table." << std::endl;
      break;
    case Opacity::userfunc:
      idfx::cout << "from a user-defined function."
                     << std::endl;
      if(!data->radiation[0]->xiFunc) {
        IDEFIX_ERROR("No opacity function has been enrolled for xi");
      }
      break;
  }
  idfx::cout << "Source term solver is ";
  switch(source_solver) {
    case ImplicitSolver::full_implicit:
      idfx::cout << "full_implicit." << std::endl;
      break;
    case ImplicitSolver::fixed_point_rad:
      idfx::cout << "fixed_point_rad." << std::endl;
      break;
    case ImplicitSolver::fixed_point_gas:
      idfx::cout << "fixed_point_gas (to test)." << std::endl;
      break;
  }

 #else
    IDEFIX_ERROR("Isothermal EOS is not compatible with radiative transfer.");
 #endif

}

void RadSource::AddRadSource(const real dt) {
  idfx::pushRegion("RadSource::AddRadSource");

  if(haveRelativistCorrection) RadSource::RelativistCorrection(dt);

  switch(source_solver) {
    case ImplicitSolver::full_implicit:
      RadSource::SourceFullImplicit(dt);
      break;
    case ImplicitSolver::fixed_point_rad:
      RadSource::SourceFixedPointRad(dt);
      break;
    case ImplicitSolver::fixed_point_gas:
      RadSource::SourceFixedPointGas(dt);
      break;
  }

  idfx::popRegion();
}



void RadSource::IrrFlux(IdefixArray3D<real> divFin) {
  idfx::pushRegion("RadSource::IrrFlux");

  auto VcGas = this->VcGas;
  IdefixArray3D<real>  dV = this->data->dV;
  IdefixArray3D<real>  A1 = this->data->A[IDIR];
  IdefixArray1D<real>  x1l = this->data->xl[IDIR];
  auto units=idfx::units;
  auto irr_type = this->irr_type;
  auto irr1D = this->irr_1D;
  IdefixArray3D<real> divFlux = divFin;
  IdefixArray3D<real> kapparho = this->kapparhoArr;
  IdefixArray3D<real> kappapArr = this->kappapArr;
  IdefixArray3D<real> tau("tau",this->data->np_tot[KDIR],this->data->np_tot[JDIR],this->data->np_tot[IDIR]);
  real kappa_irr = this->kappa_irr;

  real flux_pre = std::pow(rs/units.GetLength(),2.)*units.sigma_sb*std::pow(Ts,4.)/units.GetLength();

  if (irr_type==IradiationFlux::constant){
    column_rho->ComputeColumn(this->VcGas,RHO);
    tau = column_rho->GetColumn();
  } else if (irr_type==IradiationFlux::usertable){
    column_rho->ComputeColumn(this->VcGas,RHO);
    tau = column_rho->GetColumn();
  } else if (irr_type==IradiationFlux::userfunc){
    idefix_for("RadSourceInitKapparho",
    data->beg[KDIR], data->end[KDIR],
    data->beg[JDIR], data->end[JDIR],
    data->beg[IDIR], data->end[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
                kapparho(k,j,i) = kappapArr(k,j,i)*units.GetDensity()*units.GetLength()*VcGas(RHO,k,j,i);
    });
    column_rho->ComputeColumn(kapparho);
    tau = column_rho->GetColumn();
  }

  idefix_for("RadSourceIrrFlux",
  data->beg[KDIR], data->end[KDIR],
  data->beg[JDIR], data->end[JDIR],
  data->beg[IDIR], data->end[IDIR],
  KOKKOS_LAMBDA (int k, int j, int i) {

              real Fip,Fim;

              // Constant kappa
              if(irr_type==IradiationFlux::constant) {
                real kirr = kappa_irr*units.GetDensity()*units.GetLength();
                Fim = std::exp(-kirr*tau(k,j,i-1))*A1(k,j,i)/std::pow(x1l(i),2.);
                Fip = std::exp(-kirr*tau(k,j,i))*A1(k,j,i+1)/std::pow(x1l(i+1),2.);

              // Usertable kappa
              } else if (irr_type==IradiationFlux::usertable) {
                real logtaum = std::log10(FMAX(tau(k,j,i-1)*units.GetDensity()*units.GetLength(),1.e-15));
                Fim = pow(10.,irr1D.Get(&logtaum))*A1(k,j,i)/std::pow(x1l(i),2.);
                real logtaup = std::log10(FMAX(tau(k,j,i)*units.GetDensity()*units.GetLength(),1.e-15));
                Fip = pow(10.,irr1D.Get(&logtaup))*A1(k,j,i+1)/std::pow(x1l(i+1),2.);

              // Userfunc kappa
              } else if (irr_type==IradiationFlux::userfunc){
                Fim = std::exp(-tau(k,j,i-1))*A1(k,j,i)/std::pow(x1l(i),2.);
                Fip = std::exp(-tau(k,j,i))*A1(k,j,i+1)/std::pow(x1l(i+1),2.);
              }

              divFlux(k,j,i) = flux_pre*(Fip-Fim)/dV(k,j,i);
    });


  idfx::popRegion();

}


  void RadSource::Qviscous(IdefixArray3D<real> Qviscin) {
  idfx::pushRegion("RadSource::Qvisc");

  auto VcGas = this->VcGas;
  IdefixArray3D<real>  dV = this->data->dV;
  IdefixArray3D<real>  A1 = this->data->A[IDIR];
  IdefixArray1D<real>  x1l = this->data->xl[IDIR];
  IdefixArray1D<real>  x1 = this->data->x[IDIR];
  auto units=idfx::units;

  IdefixArray3D<real> eta1 = this->data->hydro->viscosity->eta1Arr;
  IdefixArray3D<real> eta2 = this->data->hydro->viscosity->eta2Arr;
  this->data->hydro->viscosity->viscousDiffusivityFunc(*(this->data), this->data->t, eta1, eta2);

  IdefixArray3D<real> Qvisc = this->Qvisc;

  idefix_for("RadSourceQvisc",
  data->beg[KDIR], data->end[KDIR],
  data->beg[JDIR], data->end[JDIR],
  data->beg[IDIR], data->end[IDIR],
  KOKKOS_LAMBDA (int k, int j, int i) {

              real Omegam,Omegap,rdOmegadr,mu;

              Omegam = VcGas(VX3,k,j,i)/x1(i);
              Omegap = VcGas(VX3,k,j,i+1)/x1(i+1);
              rdOmegadr = x1(i)*(Omegap-Omegam)/(x1(i+1)-x1(i));
              mu = eta1(k,j,i);

              Qvisc(k,j,i) = mu*rdOmegadr*rdOmegadr;
              //std::printf("Qvisc=%e\n",Qvisc(k,j,i));

    });


  idfx::popRegion();


}
