// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file test/Radiation/sublimation_3D_MRI/setup.cpp
//
// Last modified : 02/2026
//
// Copyright(C) by :
// - Nicolas Scepi <nicolas.scepi@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "idefix.hpp"
#include "setup.hpp"
#include "analysis.hpp"
#include "constrainedTransport.hpp"

real epsilonGlob;
real epsilonTopGlob;
real betaGlob;
real HidealGlob;
real AmMidGlob;
real LHMidGlob;
real gammaGlob;
real densityFloorGlob;
real trSmoothingGlob;
real diffCap;
Analysis *analysis;
bool haveHall;
real alphaGlob;
real alphaBetaGlob;
real muGlob;
real GGlob;
real MGlob;
real R0Glob;
real rho0Glob;
real TsubGlob;
real rhosGlob;
real rhosindexGlob;
real TwidthGlob;
real f0Glob;
real kappastarGlob;
real kappagasGlob;
real rsGlob;
real TsGlob;
real ToutGlob;

Column *columnrGlob;
Column *columnthupGlob;
Column *columnthdownGlob;

IdefixArray2D<real> *rhoInit;
std::vector<real> taus;

ConstrainedTransport<DefaultPhysics> *emf;

/*********************************************/
/**
Customized random number generator
Allow one to have consistant random numbers
generators on different architectures.
**/
/*********************************************/
real randm(void) {
    const int a    =    16807;
    const int m =    2147483647;
    static int in0 = 13763 + 2417*idfx::prank;
    int q;

    /* find random number  */
    q= (int) fmod((double) a * in0, m);
    in0=q;

    return((real) ((double) q/(double)m));
}

KOKKOS_INLINE_FUNCTION real computeDensityFloor(real R, real z, real d_floor_0, real Rin, real c0){

  real  D_return ;
  if (R>Rin){
    D_return = d_floor_0*Rin*sqrt(Rin) / (R*sqrt(R)) * Rin*Rin/(z*z+1.2*(c0*R)*(c0*R));
  }
  else{
    D_return = d_floor_0 * Rin*Rin/(z*z+1.2*(c0*Rin)*(c0*Rin));
  }
  return D_return;
}

void MyKappa(DataBlock &data, IdefixArray3D<real> &kappap, IdefixArray3D<real> &kappar) {
  IdefixArray4D<real> Vc=data.hydro->Vc;
  auto units = idfx::units;

  IdefixArray3D<real> tau;
  IdefixArray3D<real> kappa=(&data)->radiation[0]->radsource->kappapArr;
  IdefixArray3D<real> kapparho("rho",data.np_tot[KDIR],data.np_tot[JDIR],data.np_tot[IDIR]);
  idefix_for("init rho",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
    KOKKOS_LAMBDA(int k, int j, int i) {
      kapparho(k,j,i) = Vc(RHO,k,j,i)*kappa(k,j,i)*units.GetDensity()*units.GetLength();
    });

  columnrGlob->ComputeColumn(kapparho);
  tau = columnrGlob->GetColumn();

  IdefixArray1D<real> dr = data.dx[IDIR];

  real Tsub = TsubGlob;
  real rhos = rhosGlob;
  real rhos_index = rhosindexGlob;
  real Twidth = TwidthGlob;
  real f0 = f0Glob;
  real kappa_star = kappastarGlob;
  real kappa_gas = kappagasGlob;
  real mu = muGlob;

  idefix_for("MyKappa",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {

                real Tsublim = Tsub*std::pow(Vc(RHO,k,j,i)*units.GetDensity()/rhos,rhos_index);
                real T = Vc(PRS,k,j,i)/Vc(RHO,k,j,i)*units.GetKelvin()*mu;
                real f_delta = 0.2/(Vc(RHO,k,j,i)*units.GetDensity()*kappa_star*dr(i)*units.GetLength())-kappa_gas/kappa_star;
                real f_gtod;
                if ((T < Tsublim) || (tau(k,j,i) > 3.)){
                  f_gtod = f0;
                } else {
                  f_gtod = f_delta*0.25*(1.-std::tanh(std::pow((T-Tsublim)/Twidth,3.)));
                }
                f_gtod *= 1.-std::tanh(2./3.-tau(k,j,i));
                //f_gtod = 1.e-3;
                kappap(k,j,i) = kappa_star*FMAX(1.e-10,f_gtod)+kappa_gas;
                //kappar(k,j,i) = kappa_star*f_gtod+kappa_gas;
                kappar(k,j,i) = 0.;
              });
}

void Ambipolar(DataBlock& data, real t, IdefixArray3D<real> &xAin) {
  IdefixArray3D<real> xA = xAin;
  IdefixArray1D<real> x1=data.x[IDIR];
  IdefixArray1D<real> x2=data.x[JDIR];
  IdefixArray4D<real> Vc=data.hydro->Vc;
  real GM = GGlob*MGlob;
  real R0 = R0Glob;
  real epsilon = epsilonGlob;
  real etamax = diffCap*epsilon*epsilon*std::sqrt(GM/(R0*R0*R0)); // Corresponds to Rm=0.1

  auto units = idfx::units;
  real mu = muGlob;

  real xi_rad = 1.e-19; // Radioactive decay
  real K_a = 1.e-7;

  columnrGlob->ComputeColumn(data.hydro->Vc,RHO);
  columnthupGlob->ComputeColumn(data.hydro->Vc,RHO);
  columnthdownGlob->ComputeColumn(data.hydro->Vc,RHO);

  IdefixArray3D<real> Sigmar, Sigmathup, Sigmathdown;
  Sigmar = columnrGlob->GetColumn();
  Sigmathup = columnthupGlob->GetColumn();
  Sigmathdown = columnthdownGlob->GetColumn();


  idefix_for("Ambipolar",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                    real r=x1(i);
                    real th=x2(j);
                    real z=r*cos(th);
                    real R=FMAX(FABS(r*sin(th)),R0);
                    real Omega = std::sqrt(GM/std::pow(R,3.));

                    real T = Vc(PRS,k,j,i)/Vc(RHO,k,j,i)*units.GetKelvin()*mu;

                    real n_n = Vc(RHO,k,j,i)*units.GetDensity()/units.m_n;

                    // Cosmic-ray ionisation ###
                    real xi_cr = 1.e-16*(std::exp(-Sigmathup(k,j,i)*units.GetDensity()*units.GetLength()/96.)+std::exp(-Sigmathdown(k,j,i)*units.GetDensity()*units.GetLength()/96.));

                    // X-ray ionisation
                    real xi_X = std::pow(R,-2.2)*(4.0e-11*std::exp(-std::pow(Sigmar(k,j,i)*units.GetDensity()*units.GetLength()/units.m_p/3e21,0.5))+2e-14*std::exp(-pow(Sigmar(k,j,i)*units.GetDensity()*units.GetLength()/units.m_p/1.e24,0.7)));

                    // Thermal ionisation (eq. 1 of Fromang, Terquem & Blabus 2002) ###
                    real xe_th = FMIN(1.,6.47e-13*std::sqrt(K_a/1.e-7)*std::pow((T/1.e3),0.75)*std::sqrt(2.4e15/n_n)*std::exp(-25188/T)/1.15e-11);

                    // Ionisation due to UV flux
                    real xe_UV = 2.e-5*std::exp(-(Sigmar(k,j,i)*units.GetDensity()*units.GetLength()/0.03));

                    real xi = xi_rad+xi_cr+xi_X;

                    real alpha_dr = 3.e-6/std::sqrt(T);

                    real xe_nth = std::sqrt(xi/(n_n*alpha_dr));

                    real xe = FMIN(1.,xe_th+xe_nth+xe_UV);

                    real gammai = 1.3e-9/(units.m_n+units.m_p);

                    real B2 = (Vc(BX1,k,j,i)*Vc(BX1,k,j,i)+Vc(BX2,k,j,i)*Vc(BX2,k,j,i)+Vc(BX3,k,j,i)*Vc(BX3,k,j,i))*units.GetVelocity()*units.GetVelocity()*units.GetDensity();

                    // Eq .14c of Lesur, Kunz & Fromang 2014
                    real eta = B2/(gammai*pow(Vc(RHO,k,j,i)*units.GetDensity(),2.)*xe);

                    // Use a diffusivity cap that leads to a constant CFL in radius
                        real etamax_loc = etamax*x1(i)*x1(i)*units.GetVelocity()*units.GetLength();

                    if(eta>etamax_loc) xA(k,j,i) = etamax_loc/B2*(units.GetDensity()/units.GetTime());
                    else xA(k,j,i) = eta/B2*(units.GetDensity()/units.GetTime());
                    //if (i==2 && j==2 && k==3) std::printf("xA=%e, eta=%e, rho=%e, T=%e, xi_cr=%e, xi_X=%e, xe_th=%e, xe_UV=%e, alpha_dr=%e, Sigmar=%e\n",xA(k,j,i),eta,Vc(RHO,k,j,i)*units.GetDensity(),T,xi_cr,xi_X,xe_th,xe_UV,alpha_dr,Sigmar(k,j,i));

              });



}

void Resistivity(DataBlock& data, real t, IdefixArray3D<real> &etain) {
  IdefixArray3D<real> eta = etain;
  IdefixArray1D<real> x1=data.x[IDIR];
  IdefixArray1D<real> x2=data.x[JDIR];
  IdefixArray4D<real> Vc=data.hydro->Vc;

  auto units = idfx::units;
  real epsilon = epsilonGlob;
  real mu = muGlob;
  real R0 = R0Glob;
  real GM = GGlob*MGlob;
  real etamax = diffCap*epsilon*epsilon*std::sqrt(GM/(R0*R0*R0)); // Corresponds to Rm=0.1

  real xi_rad = 1.e-19; // Radioactive decay
  real K_a = 1.e-7;

  columnrGlob->ComputeColumn(data.hydro->Vc,RHO);
  columnthupGlob->ComputeColumn(data.hydro->Vc,RHO);
  columnthdownGlob->ComputeColumn(data.hydro->Vc,RHO);

  IdefixArray3D<real> Sigmar, Sigmathup, Sigmathdown;
  Sigmar = columnrGlob->GetColumn();
  Sigmathup = columnthupGlob->GetColumn();
  Sigmathdown = columnthdownGlob->GetColumn();

  idefix_for("Resistivity",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                    real r=x1(i);
                    real th=x2(j);
                    real z=r*cos(th);
                    real R=r*sin(th);

                    real Omega = std::sqrt(GM/std::pow(R,3.));
                    real cs = std::sqrt(Vc(PRS,k,j,i)/Vc(RHO,k,j,i));

                    real T = Vc(PRS,k,j,i)/Vc(RHO,k,j,i)*units.GetKelvin()*mu;

                    real n_n = Vc(RHO,k,j,i)*units.GetDensity()/units.m_n;

                    // Cosmic-ray ionisation ###
                    real xi_cr = 1.e-16*(std::exp(-Sigmathup(k,j,i)*units.GetDensity()*units.GetLength()/96.)+std::exp(-Sigmathdown(k,j,i)*units.GetDensity()*units.GetLength()/96.));

                    // X-ray ionisation
                    real xi_X = std::pow(R/R0,-2.2)*(4.0e-11*std::exp(-std::pow(Sigmar(k,j,i)*units.GetDensity()*units.GetLength()/units.m_p/3e21,0.5))+2e-14*std::exp(-std::pow(Sigmar(k,j,i)*units.GetDensity()*units.GetLength()/units.m_p/1.e24,0.7)));

                    // Thermal ionisation (eq. 1 of Fromang, Terquem & Blabus 2002) ###
                    real xe_th = FMIN(1.,6.47e-13*std::sqrt(K_a/1.e-7)*std::pow((T/1.e3),0.75)*std::sqrt(2.4e15/n_n)*std::exp(-25188./T)/1.15e-11);
                    //real xe_th = FMIN(1.e-3,6.47e-13*std::sqrt(K_a/1.e-7)*std::pow((T/1.e3),0.75)*std::sqrt(2.4e15/n_n)*std::exp(-25188./T)/1.15e-11*std::exp(-z*z/(2.*h*h)));
                    //real xe_th = 1.e-5;

                    // Ionisation due to UV flux
                    real xe_UV = 2.e-5*std::exp(-(Sigmar(k,j,i)*units.GetDensity()*units.GetLength()/0.03));

                    real xi = xi_rad+xi_cr+xi_X;

                    real alpha_dr = 3.e-6/std::sqrt(T);

                    real xe_nth = std::sqrt(xi/(n_n*alpha_dr));

                    real xe = FMIN(1.,xe_th+xe_nth+xe_UV);
                    //real xe = 1.e-5;

                    real sigmave = 8.28e-9*std::sqrt(T/100.);

                    // Eq .14a of Lesur, Kunz & Fromang 2014
                    //eta(k,j,i) = units.c*units.c*units.m_e*sigmave/(4.*M_PI*units.e*units.e*xe);
                    eta(k,j,i) = units.c*units.c*units.m_e*sigmave/(4.*M_PI*units.e*units.e*xe)*(units.GetTime()/(units.GetLength()*units.GetLength()));

                    // Use a diffusivity cap that leads to a constant CFL in radius
                        real etamax_loc = etamax*x1(i)*x1(i)*units.GetVelocity()*units.GetLength();
                    if(eta(k,j,i)>etamax_loc) eta(k,j,i) = etamax_loc;
                    //eta(k,j,i) = 5.e-5;

                    //real B2 = (Vc(BX1,k,j,i)*Vc(BX1,k,j,i)+Vc(BX2,k,j,i)*Vc(BX2,k,j,i)+Vc(BX3,k,j,i)*Vc(BX3,k,j,i))*units.GetVelocity()*units.GetVelocity()*units.GetDensity();

                    //if (i==2 && j==2 && k==3) std::printf("eta=%e, Rm=%e, rho=%e, T=%e, xi_cr=%e, xi_X=%e, xe_th=%e, xe_nth=%e, xe_UV=%e, alpha_dr=%e, B2=%e\n",eta(k,j,i),cs*cs/(eta(k,j,i)*Omega),Vc(RHO,k,j,i)*units.GetDensity(),T,xi_cr,xi_X,xe_th,xe_nth,xe_UV,alpha_dr,B2);
                    //if (i==10 && k==0) std::printf("eta=%e, rho=%e, T=%e, xi_cr=%e, xi_X=%e, xe_th=%e, xe_UV=%e, alpha_dr=%e\n",eta(k,j,i),Vc(RHO,k,j,i)*units.GetDensity(),T,xi_cr,xi_X,xe_th,xe_UV,alpha_dr);

              });

}

void MySourceTerm(Fluid<DefaultPhysics> *hydro, const real t, const real dtin) {
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray4D<real> Uc = hydro->Uc;
  auto data = hydro->data;
  auto x1 = data->x[IDIR];
  auto x2 = data->x[JDIR];
  IdefixArray2D<real> rhoEq = *rhoInit;

  real epsilonTop = epsilonTopGlob;
  real epsilon = epsilonGlob;

  real Omega = std::sqrt(GGlob*MGlob)*std::pow(R0Glob,-1.5);

  real tauGlob=1.0/Omega;
  real tauWind=1e-2/Omega;
  real tauInner=0.1/Omega;
  real gamma_m1=gammaGlob-1.0;
  real dt=dtin;
  real Hideal=HidealGlob;
  real R0=R0Glob;
  real trSmoothing = trSmoothingGlob;
  real alpha=alphaGlob;
  real GM = GGlob*MGlob;

  idefix_for("MySourceTerm",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real r=x1(i);
                real th=x2(j);
                real z=r*cos(th);
                real R=r*sin(th);
                real Ri = FMAX(R0,R);
                real Vk=std::sqrt(GM/Ri);

                real Zh = FABS(z/R)/epsilon;
                real Tdisk = std::pow(epsilon*Vk,2.);  // Factor of GM to take into account vk not 1 at R=1
                real Tcorona = std::pow(epsilonTop*Vk,2.);  // Factor of GM to take into account vk not 1 at R=1
                //if(x1(i) < 1.5) cscorona = csdisk;
                real f = tanh((Zh-Hideal)/trSmoothing);
                real Teff=0.5*(Tdisk+Tcorona)+0.5*(Tcorona-Tdisk)*f;

                real tau= 0.5*(tauGlob+tauWind)+0.5*(tauWind-tauGlob)*f;
                tau *= pow(Ri/R0,1.5);
                // Cooling /heatig function
                real Ptarget = Teff*Vc(RHO,k,j,i);

                Uc(ENG,k,j,i) += -dt*(Vc(PRS,k,j,i)-Ptarget)/(tau*gamma_m1);

});


}



void InternalBoundary(DataBlock& data, const real t) {
  IdefixArray4D<real> Vc = data.hydro->Vc;
  IdefixArray4D<real> dustVc;
  bool haveDust = data.haveDust;
  if(haveDust) {
   dustVc= data.dust[0]->Vc;
  }
  IdefixArray4D<real> Vs = data.hydro->Vs;
  IdefixArray1D<real> x1=data.x[IDIR];
  IdefixArray1D<real> x2=data.x[JDIR];

  auto units = idfx::units;

  real Rin=R0Glob;
  real GM = GGlob*MGlob;

  real vAmax = 2.*std::sqrt(GM/Rin);
  real densityFloor0 = densityFloorGlob;
  real epsilon=epsilonGlob;

  idefix_for("InternalBoundary",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      real R=x1(i)*sin(x2(j));
      real z=x1(i)*cos(x2(j));
      real zh = FABS(z/R)/epsilon;

      real b2=EXPAND(Vc(BX1,k,j,i)*Vc(BX1,k,j,i) , +Vc(BX2,k,j,i)*Vc(BX2,k,j,i), +Vc(BX3,k,j,i)*Vc(BX3,k,j,i) ) ;
      real va2=b2/Vc(RHO,k,j,i);
      real myMax=vAmax;
      //if(x1(i)<1.1) myMax=myMax/50.0;
      if(va2>myMax*myMax) {
        real oldrho = Vc(RHO,k,j,i);
        real T = Vc(PRS,k,j,i)/Vc(RHO,k,j,i);
        Vc(RHO,k,j,i) = b2/(myMax*myMax);
        Vc(PRS,k,j,i) = T*Vc(RHO,k,j,i);
        Vc(VX1,k,j,i) *= oldrho/Vc(RHO,k,j,i);
        Vc(VX2,k,j,i) *= oldrho/Vc(RHO,k,j,i);
        Vc(VX3,k,j,i) *= oldrho/Vc(RHO,k,j,i);
      }
      real densityFloor = FMAX(densityFloor0,computeDensityFloor(R,z,densityFloor0,Rin,epsilon))/units.GetDensity();
      if(Vc(RHO,k,j,i) < densityFloor) {
        real oldrho = Vc(RHO,k,j,i);
        real T= Vc(PRS,k,j,i)/Vc(RHO,k,j,i);
        Vc(RHO,k,j,i)=densityFloor;
        Vc(VX1,k,j,i) *= oldrho/densityFloor;
        Vc(VX2,k,j,i) *= oldrho/densityFloor;
        Vc(VX3,k,j,i) *= oldrho/densityFloor;
      }

      /*
        real R = x1(i)*sin(x2(j));
        if(R<1.0) {
            Vc(VX1,k,j,i) = ZERO_F;
            Vc(VX2,k,j,i) = ZERO_F;
            Vc(VX3,k,j,i) = R;
        }*/
    });

}
// User-defined boundaries
void UserdefBoundary(Fluid<DefaultPhysics> *hydro, int dir, BoundarySide side, real t) {

    if( (dir==IDIR) && (side == left)) {
        IdefixArray4D<real> Vc = hydro->Vc;
        IdefixArray4D<real> Vs = hydro->Vs;
        auto *data = hydro->data;
        IdefixArray1D<real> x1 = data->x[IDIR];
        IdefixArray1D<real> x2 = data->x[JDIR];

        int ighost = data->nghost[IDIR];
        real Rin = R0Glob;
        real Omega = std::sqrt(GGlob*MGlob)*std::pow(Rin,-1.5);
        real csdisk = std::sqrt(GGlob*MGlob)*epsilonGlob/sqrt(Rin);
        real cscorona = std::sqrt(GGlob*MGlob)*epsilonTopGlob/sqrt(Rin);

        hydro->boundary->BoundaryFor("UserDefX1",dir,side,
            KOKKOS_LAMBDA (int k, int j, int i) {
                real R=x1(i)*sin(x2(j));
                real z=x1(i)*cos(x2(j));

                Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost);

                //Vc(PRS,k,j,i) = Vc(RHO,k,j,i)*csdisk*csdisk;
                Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost);

                if(Vc(VX1,k,j,ighost)>=ZERO_F) Vc(VX1,k,j,i) = -Vc(VX1,k,j,2*ighost-i-1);
                       else Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost);
                Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost);

                Vc(VX3,k,j,i) = Omega*R;
                #if DIMENSIONS < 3
                Vc(BX3,k,j,i) = - Vc(BX3,k,j,2*ighost-i-1);
                #endif

            });
      hydro->boundary->BoundaryForX2s("UserDefX2s",dir,side,
        KOKKOS_LAMBDA (int k, int j, int i) {
            Vs(BX2s,k,j,i) = Vs(BX2s,k,j,ighost);
          });
      #if DIMENSIONS == 3
      hydro->boundary->BoundaryForX3s("UserDefX3s",dir,side,
        KOKKOS_LAMBDA (int k, int j, int i) {
            Vs(BX3s,k,j,i) = -Vs(BX3s,k,j,2*ighost-i-1);
          });
      #endif
    }

    if( (dir==IDIR) && (side == right)) {
        IdefixArray4D<real> Vc = hydro->Vc;
        IdefixArray4D<real> Vs = hydro->Vs;
        auto *data = hydro->data;
        IdefixArray1D<real> x1 = data->x[IDIR];
        IdefixArray1D<real> x2 = data->x[JDIR];

        int ighost = data->end[IDIR]-1;
        real Rin = R0Glob;
        real Omega = std::sqrt(GGlob*MGlob)*std::pow(Rin,-1.5);
        real csdisk = std::sqrt(GGlob*MGlob)*epsilonGlob/sqrt(Rin);
        real cscorona = std::sqrt(GGlob*MGlob)*epsilonTopGlob/sqrt(Rin);

        hydro->boundary->BoundaryFor("UserDefX1",dir,side,
            KOKKOS_LAMBDA (int k, int j, int i) {
                real R=x1(i)*sin(x2(j));
                real z=x1(i)*cos(x2(j));

                Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost);
                Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost);

                if(Vc(VX1,k,j,ighost)<=ZERO_F) Vc(VX1,k,j,i) = 0.0;
                       else Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost);
                Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost);
                //real Rmin = FMAX(0.3,R);

                //Vc(VX3,k,j,i) = 1.0/sqrt(Rmin) * sqrt( Rmin / sqrt(Rmin*Rmin + z*z));
                Vc(VX3,k,j,i) = Vc(VX3,k,j,ighost);
                #if DIMENSIONS < 3
                Vc(BX3,k,j,i) = - Vc(BX3,k,j,2*ighost-i+1);
                #endif
                //Vc(BX3,k,j,i) = Vc(BX3,k,j,ighost);

            });
      hydro->boundary->BoundaryForX2s("UserDefX2s",dir,side,
        KOKKOS_LAMBDA (int k, int j, int i) {
            Vs(BX2s,k,j,i) = Vs(BX2s,k,j,ighost);
          });
      #if DIMENSIONS == 3
      hydro->boundary->BoundaryForX3s("UserDefX3s",dir,side,
        KOKKOS_LAMBDA (int k, int j, int i) {
            Vs(BX3s,k,j,i) = 0.0;
          });
      #endif

    }
}

void UserdefBoundaryRad(Fluid<RadiationPhysics> *radiation, int dir, BoundarySide side, real t) {
  IdefixArray4D<real> Vc = radiation->Vc;
  auto *data = radiation->data;
  auto units=idfx::units;

  real Tout = ToutGlob;

  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  if(dir==IDIR) {
    int ighost,nxi,iend,ibeg;
    if(side == left) {
      ighost = data->nghost[IDIR];
      ibeg = 0;
      iend = data->beg[IDIR];
      idefix_for("UserDefBoundaryRad",
        0, data->np_tot[KDIR],
        0, data->np_tot[JDIR],
        ibeg, iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vc(ER,k,j,i) = units.ar*std::pow(Tout,4.)/units.GetEnergy();
          if (Vc(FR1,k,j,ighost) >=ZERO_F){
            Vc(FR1,k,j,i) = ZERO_F;
          } else {
            Vc(FR1,k,j,i) = Vc(FR1,k,j,ighost);
          }
          Vc(FR2,k,j,i) = Vc(FR2,k,j,ighost);
          Vc(FR3,k,j,i) = Vc(FR3,k,j,ighost);
        });
    } else if (side==right){
      ighost = data->nghost[IDIR];
      nxi = data->np_int[IDIR];
      ibeg = data->end[IDIR];
      iend =data->np_tot[IDIR];
      idefix_for("UserDefBoundaryRad",
        0, data->np_tot[KDIR],
        0, data->np_tot[JDIR],
        ibeg, iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vc(ER,k,j,i) = units.ar*std::pow(Tout,4.)/units.GetEnergy();
          if (Vc(FR1,k,j,ighost+nxi-1) <=ZERO_F){
            Vc(FR1,k,j,i) = ZERO_F;
          } else {
            Vc(FR1,k,j,i) = Vc(FR1,k,j,ighost+nxi-1);
          }
          Vc(FR2,k,j,i) = Vc(FR2,k,j,ighost+nxi-1);
          Vc(FR3,k,j,i) = Vc(FR3,k,j,ighost+nxi-1);
        });
    }
  }
}

void EmfBoundary(DataBlock& data, const real t) {
    IdefixArray3D<real> Ex1 = data.hydro->emf->ex;
    IdefixArray3D<real> Ex2 = data.hydro->emf->ey;
    IdefixArray3D<real> Ex3 = data.hydro->emf->ez;
    if(data.lbound[IDIR] == userdef) {

        int ighost = data.nghost[IDIR];

        // Do not permit poloidal field to enter the seed
        idefix_for("EMFBoundary",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,ighost+1,
                    KOKKOS_LAMBDA (int k, int j, int i) {
            Ex3(k,j,i) = ZERO_F;
            #if DIMENSIONS == 3
              Ex2(k,j,i) = ZERO_F;
            #endif
        });
    }
}

void FluxBoundary(DataBlock & data, int dir, BoundarySide side, const real t) {
    IdefixArray4D<real> Flux = data.hydro->FluxRiemann[dir];
    if( dir==IDIR && side == left) {
        int iref = data.beg[IDIR];

        idefix_for("FluxBoundLeft",
                    data.beg[KDIR], data.end[KDIR],
                    data.beg[JDIR], data.end[JDIR],
          KOKKOS_LAMBDA (int k, int j) {
            if(Flux(RHO, k, j, iref) > 0.0) {
              Flux(RHO, k, j, iref) = 0.0; // Cancel incoming mass flux.
            }
          });
    }

}

void CoarsenFunction(DataBlock &data) {
  IdefixArray2D<int> coarseningLevel = data.coarseningLevel[KDIR];
  IdefixArray1D<real> th = data.x[JDIR];
  idefix_for("set_coarsening", 0, data.np_tot[JDIR], 0, data.np_tot[IDIR],
      KOKKOS_LAMBDA(int j,int i) {
        int c = 1.0/fabs(sin(th(j)));
        if(c>5) c = 5;
        coarseningLevel(j,i) = c;

      });
}


void ComputeUserVars(DataBlock & data, UserDefVariablesContainer &variables) {

  // custom scratch array
  IdefixArray3D<real> scrh("Scratch", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);
  IdefixArray3D<real> scrh2("Scratch", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);
  IdefixArray3D<real> xH;
  IdefixArray3D<real> Sigmar,Sigmathup,Sigmathdown;
  IdefixArray3D<real> kappap, kappar;

  Sigmar = columnrGlob->GetColumn();
  Sigmathup = columnthupGlob->GetColumn();
  Sigmathdown = columnthdownGlob->GetColumn();

  Kokkos::deep_copy(variables["columnr"],Sigmar);
  Kokkos::deep_copy(variables["columnthup"],Sigmathup);
  Kokkos::deep_copy(variables["columnthdown"],Sigmathdown);
  Kokkos::deep_copy(variables["kappap"],data.radiation[0]->radsource->kappapArr);
  Kokkos::deep_copy(variables["kappar"],data.radiation[0]->radsource->kapparArr);

  // Ask for a computation of xA ambipolar in this scratch array
  Ambipolar(data, data.t, scrh);
  Resistivity(data, data.t, scrh2);

  // Mirror data on Host
  DataBlockHost d(data);

  // Sync it
  d.SyncFromDevice();

  // Make references to the user-defined arrays (variables is a container of IdefixHostArray3D)
  // Note that the labels should match the variable names in the input file
  IdefixHostArray3D<real> xA  = variables["xA"];
  IdefixHostArray3D<real> etaO  = variables["etaO"];

  Kokkos::deep_copy(xA,scrh);
  Kokkos::deep_copy(etaO,scrh2);

  IdefixHostArray1D<real> x1=d.x[IDIR];
  IdefixHostArray1D<real> x2=d.x[JDIR];
  IdefixHostArray4D<real> Vc=d.Vc;

  for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
    for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
      for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
        variables["InvDt"](k,j,i) = d.InvDt(k,j,i);
      }
    }
  }
}

void analysisFunction(DataBlock& data) {
  analysis->PerformAnalysis(data);
}


void ComputeRho(DataBlock &data) {
  real Rin=R0Glob;
  real alpha = alphaGlob; // Density power law
  real GM = GGlob*MGlob;
  auto units = idfx::units;

  DataBlockHost d(data);
  // Compute vertical equilibrium following Bai & Stone (2017) (BS17)
  GridHost gh(*data.mygrid);
  gh.SyncFromDevice();
  IdefixHostArray1D<real> logf = IdefixHostArray1D<real>("fvert", gh.np_tot[JDIR]);
  int jmid = gh.np_tot[JDIR]/2;
  logf(jmid) = 0;

  for(int j = jmid+1 ; j < gh.np_tot[JDIR] ; j++) {
    real th = gh.xl[JDIR](j);
    real dth = 0.5*(gh.dx[JDIR](j-1)+gh.dx[JDIR](j));
    real Zh = 1/FABS(tan(th))/epsilonGlob;
    real Vk=std::sqrt(GM)/pow(Rin,0.5);
    real Tdisk = std::pow(epsilonGlob*Vk,2.);  // Factor of GM to take into account vk not 1 at R=1
    real Tcorona = std::pow(epsilonTopGlob*Vk,2.);
    //if(x1(i) < 1.5) cscorona = csdisk;
    real Teff=0.5*(Tdisk+Tcorona)+0.5*(Tcorona-Tdisk)*tanh((Zh-HidealGlob)/trSmoothingGlob);
    // BS17 eq. 14
    logf(j) = logf(j-1)+dth*(cos(th)/sin(th)*GM/Teff/Rin-(alpha+1));
    //idfx::cout << "j=" << j << " th-pi/2=" << th - M_PI/2 << " Zh="<< Zh << " Teff=" << Teff << " logf=" << logf(j) << " log(sin(th))/g=" << log(sin(th))/Teff << std::endl;
  }

  // That's for the other side, by symmetry accross the midplane
  for(int j = 0 ; j < jmid ; j++) {
    logf(j) = logf(2*jmid-j-1);
  }

  rhoInit = new IdefixArray2D<real>("RhoInit",d.np_tot[JDIR],d.np_tot[IDIR]);
  IdefixHostArray2D<real> rhoH("RhoHost",d.np_tot[JDIR],d.np_tot[IDIR]);

  for(int j = 0; j < d.np_tot[JDIR] ; j++) {
    int jglob = j + d.gbeg[JDIR] - d.beg[JDIR];
    for(int i = 0; i < d.np_tot[IDIR] ; i++) {
      real r=d.x[IDIR](i);
      real th=d.x[JDIR](j);
      real z=r*cos(th);
      real R=r*sin(th);
      real Ri = FMAX(R,Rin);
      real Zh = FABS(z/R)/epsilonGlob;

      real Vk=std::sqrt(GM/Rin);
      real Tdisk = std::pow(epsilonGlob*Vk,2.);  // Factor of GM to take into account vk not 1 at R=1
      real Tcorona = std::pow(epsilonTopGlob*Vk,2.);
      //if(x1(i) < 1.5) cscorona = csdisk;
      real Teff=0.5*(Tdisk+Tcorona)+0.5*(Tcorona-Tdisk)*tanh((Zh-HidealGlob)/trSmoothingGlob);

      // BS17 F(theta)=f(theta)*g(theta)
      // here logF is log(F) and g is Teff
      // We normalise f so that it is=1 in the disc midaplne
      real f = exp(logf(jglob)) *epsilonGlob*epsilonGlob*Vk*Vk / Teff;


      rhoH(j,i) = pow(Ri/R0Glob,-alpha) * f * rho0Glob/units.GetDensity();

      real densityFloor = FMAX(densityFloorGlob,computeDensityFloor(R,z,densityFloorGlob,Rin,epsilonGlob))/units.GetDensity();
      if(rhoH(j,i) < densityFloor) {
        rhoH(j,i) = densityFloor;
      }
    }}
  Kokkos::deep_copy(*rhoInit,rhoH);
}

// Default constructor


// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  // Copy the pointer to the emf function for future usage
  emf = data.hydro->emf.get();
  // Set the function for userdefboundary
  data.hydro->EnrollUserDefBoundary(&UserdefBoundary);

  data.hydro->EnrollAmbipolarDiffusivity(&Ambipolar);

  if(data.haveGridCoarsening) {
    data.EnrollGridCoarseningLevels(&CoarsenFunction);
  }

  if(data.hydro->resistivityStatus.status != HydroModuleStatus::Disabled) {
    data.hydro->EnrollOhmicDiffusivity(&Resistivity);
  }
  data.hydro->EnrollInternalBoundary(&InternalBoundary);
  //data.hydro->EnrollEmfBoundary(&EmfBoundary);
  //data.hydro->EnrollFluxBoundary(&FluxBoundary);
  //data.hydro->EnrollUserSourceTerm(&MySourceTerm);

  gammaGlob=data.hydro->eos->GetGamma();
  alphaGlob = input.GetOrSet<real>("Setup","alpha",0,1.5);  // density power law
  alphaBetaGlob = input.GetOrSet<real>("Setup","alphaBeta",0,0.0); // beta power law
  epsilonGlob = input.Get<real>("Setup","epsilon",0);
  epsilonTopGlob = input.Get<real>("Setup","epsilonTop",0);
  betaGlob = input.Get<real>("Setup","beta",0);
  HidealGlob = input.Get<real>("Setup","Hideal",0);
  diffCap = input.Get<real>("Setup","diffusivityCap",0);
  if(haveHall) LHMidGlob = input.GetOrSet<real>("Setup","LH",0,0.0);
  rho0Glob = input.Get<real>("Setup","rho0",0);

  densityFloorGlob = input.Get<real>("Setup","density_floor",0);
  trSmoothingGlob = input.Get<real>("Setup","transitionSmoothing",0);

  muGlob = input.Get<real>("Hydro","mu",0);

  GGlob = input.Get<real>("Gravity","gravCst",0);
  MGlob = input.Get<real>("Gravity","Mcentral",0);
  R0Glob = input.Get<real>("Setup","R0",0);

  columnrGlob = new Column(IDIR,1,&data);
  columnthupGlob = new Column(JDIR,1,&data);
  columnthdownGlob = new Column(JDIR,-1,&data);

  output.EnrollUserDefVariables(&ComputeUserVars);
  // assume disc surface at 6.5 h
  analysis = new Analysis(grid, data,std::string("profile.dat"), HidealGlob*epsilonGlob);
  output.EnrollAnalysis(&analysisFunction);

  rsGlob=input.Get<real>("Rad","irr",1);
  TsGlob=input.Get<real>("Rad","irr",2);
  TsubGlob = input.Get<real>("Setup","Tsub",0);
  rhosGlob = input.Get<real>("Setup","rhos",0);
  rhosindexGlob = input.Get<real>("Setup","rhos_index",0);
  TwidthGlob = input.Get<real>("Setup","Twidth",0);
  f0Glob = input.Get<real>("Setup","f0",0);
  kappastarGlob = input.Get<real>("Setup","kappa_star",0);
  kappagasGlob = input.Get<real>("Setup","kappa_gas",0);

  densityFloorGlob = input.Get<real>("Setup","density_floor",0);
  data.radiation[0]->EnrollKappa(&MyKappa);
  data.radiation[0]->EnrollUserDefBoundary(&UserdefBoundaryRad);

  ComputeRho(data);

}



Setup::~Setup() {
  delete rhoInit;
  delete analysis;
}
// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);
    auto units = idfx::units;

    // Make vector potential
    IdefixHostArray4D<real> A = IdefixHostArray4D<real>("Setup_VectorPotential", 3, data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);

    // Get the equilibrium density profile
    IdefixHostArray2D<real> rhoH("RhoHost",d.np_tot[JDIR],d.np_tot[IDIR]);
    Kokkos::deep_copy(rhoH,*rhoInit);

    real Rin=R0Glob;
    real alpha = alphaGlob; // MINUS Density power law
    real gammaB =0.5*(-alpha-1-alphaBetaGlob); // MAgnetic field power law

    real GM = GGlob*MGlob;

    real B0 = epsilonGlob*sqrt(2.0/betaGlob)*std::sqrt(rho0Glob/units.GetDensity())*std::sqrt(GM/Rin);

    int nSpecies = data.dust.size();

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                real r=d.x[IDIR](i);
                real th=d.x[JDIR](j);
                real z=r*cos(th);
                real R=r*sin(th);
                real Ri = FMAX(R,Rin);
                real Zh = FABS(z/R)/epsilonGlob;


                real Vk=std::sqrt(GM/Rin);
                real Tdisk = std::pow(epsilonGlob*Vk,2.);  // Factor of GM to take into account vk not 1 at R=1
                real Tcorona = std::pow(epsilonTopGlob*Vk,2.);
                //if(x1(i) < 1.5) cscorona = csdisk;
                real Teff=0.5*(Tdisk+Tcorona)+0.5*(Tcorona-Tdisk)*tanh((Zh-HidealGlob)/trSmoothingGlob);

                d.Vc(RHO,k,j,i) = rhoH(j,i);

                d.Vc(PRS,k,j,i) = Teff*Rin/Ri * d.Vc(RHO,k,j,i);
                real T = d.Vc(PRS,k,j,i)/d.Vc(RHO,k,j,i)*units.GetKelvin();

                // BS17 eq. 15
                d.Vc(VX3,k,j,i) = std::sqrt(GM)*sqrt((1-(alpha+1)*Teff)/Ri);

                // Enforce solid body rotation above the "seed"
                if(R<Rin) d.Vc(VX3,k,j,i) = std::sqrt(GM)*sqrt((1-(alpha+1)*Teff))*R/Rin;

                d.Vc(VX1,k,j,i) = d.Vc(VX3,k,j,i)*1e-1*(0.5-randm());
                d.Vc(VX2,k,j,i) = ZERO_F;

                // Vector potential on the corner
                real s=sin(d.xl[JDIR](j));
                R=d.xl[IDIR](i) * s;

                A(IDIR,k,j,i) = ZERO_F;
                A(JDIR,k,j,i) = ZERO_F;

                #ifdef EVOLVE_VECTOR_POTENTIAL
                  if(R>Rin) {
                    d.Ve(AX3e,k,j,i) = B0*(
                                      1/(gammaB+2)*(pow(R,gammaB+1) - pow(Rin,gammaB+2)/R)
                                      + Rin*Rin/(2.0*R));
            //d.Ve(AX3e,k,j,i) = B0*(pow(Rin,m+2.0)/R * (-1.0/(m+2.0)) + pow(R,m+1.0)/(m+2.0));
                  }
                  else {
                    d.Ve(AX3e,k,j,i) = B0*R/2.0;
            //d.Ve(AX3e,k,j,i) = 0.0;
                  }
                #else
                  if(R>Rin) {
                    A(KDIR,k,j,i) = B0*(
                                      Rin/(gammaB+2)*(pow(R/Rin,gammaB+1) - Rin/R)
                                      + Rin*Rin/(2.0*R));
                  }
                  else {
                    A(KDIR,k,j,i) = B0*R/2.0;
                  }
                #endif

            }
        }
    }

    // Make the field from the vector potential
    #ifndef EVOLVE_VECTOR_POTENTIAL
      d.MakeVsFromAmag(A);
    #endif


    // Send it all, if needed
    d.SyncToDevice();
}
