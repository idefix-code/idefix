// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file test/Radiation/novae_like_3D/setup.cpp
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
    D_return = d_floor_0 / (R*sqrt(R)) * 1.0/(z*z+1.2*(c0*R)*(c0*R));
  }
  else{
    D_return = d_floor_0 / (Rin*sqrt(Rin)) * 1.0/(z*z+1.2*(c0*Rin)*(c0*Rin));
  }
  if (D_return < 1.0e-9){
    D_return = 1e-9;
  }
  return D_return;
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
  real tauGlob=1.0;
  real tauWind=1e-2;
  real tauInner=0.1;
  real gamma_m1=gammaGlob-1.0;
  real dt=dtin;
  real Hideal=HidealGlob;
  //real tauVel=0.5;
  real R0=1.0;
  real trSmoothing = trSmoothingGlob;
  real alpha=alphaGlob;

  idefix_for("MySourceTerm",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real r=x1(i);
                real th=x2(j);
                real z=r*cos(th);
                real R=r*sin(th);
                real Ri = FMAX(R0,R);
                real Vk=1.0/pow(Ri,0.5);

                real Zh = FABS(z/R)/epsilon;
                real Tdisk = epsilon*epsilon/Ri;
                real Tcorona = epsilonTop*epsilonTop/Ri;
                //if(x1(i) < 1.5) cscorona = csdisk;
                real f = tanh((Zh-Hideal)/trSmoothing);
                real Teff=0.5*(Tdisk+Tcorona)+0.5*(Tcorona-Tdisk)*f;

                real tau= 0.5*(tauGlob+tauWind)+0.5*(tauWind-tauGlob)*f;
                tau *= pow(Ri,1.5);
                // Cooling /heatig function
                real Ptarget = Teff*Vc(RHO,k,j,i);

                Uc(ENG,k,j,i) += -dt*(Vc(PRS,k,j,i)-Ptarget)/(tau*gamma_m1);

                // Spatial Damping function
                real lambda = 1/tauInner * (FMAX((1.2*R0-r)/0.2,0.0));
                real rhoTarget = rhoEq(j,i);

                real vx3Target = sqrt((1-(alpha+1)*Teff)/Ri);

                // Enforce solid body rotation above the "seed"
                if(R<R0) vx3Target = sqrt((1-(alpha+1)*Teff))*R/R0;

                // With fargo, there is no mean vx3!
                //vx3Target = 0.0;

                // relaxation on all components
                real rho = Uc(RHO,k,j,i);
                real vx1 = Uc(MX1,k,j,i)/rho;
                real vx2 = Uc(MX2,k,j,i)/rho;
                real vx3 = Uc(MX3,k,j,i)/rho;
                real bx1 = Uc(BX1,k,j,i);
                real bx2 = Uc(BX2,k,j,i);
                real bx3 = Uc(BX3,k,j,i);

                real ek = 0.5*rho*(vx1*vx1+vx2*vx2+vx3*vx3);
                real em = 0.5*(bx1*bx1+bx2*bx2+bx3*bx3);

                real prs = gamma_m1*(Uc(ENG,k,j,i) - ek - em);
                real T = prs/rho;


                rho -= lambda*(Vc(RHO,k,j,i)-rhoTarget)*dt;
                vx1 -= lambda*Vc(VX1,k,j,i)*dt;
                vx2 -= lambda*Vc(VX2,k,j,i)*dt;
                vx3 -= lambda*(Vc(VX3,k,j,i)-vx3Target)*dt;

                prs = T*rho;
                ek = 0.5*rho*(vx1*vx1+vx2*vx2+vx3*vx3);

                Uc(RHO,k,j,i) = rho;
                Uc(MX1,k,j,i) = rho*vx1;
                Uc(MX2,k,j,i) = rho*vx2;
                Uc(MX3,k,j,i) = rho*vx3;
                Uc(ENG,k,j,i) = ek+em+prs/gamma_m1;

                // inner shell relaxation
                /*
                if(R<Rin) {
                  real rhoTarget = 1.0/(R0*sqrt(R0))  * exp(1.0/ Tdisk * (1.0/sqrt(R0*R0+z*z)-1.0/R0));
                  real densityFloor = computeDensityFloor(R,z,densityFloor0,Rin,epsilon);
                  if(rhoTarget < densityFloor) rhoTarget = densityFloor;

                  real vx3Target = 1.0/sqrt(R0) * sqrt( FMAX(R0 / sqrt(R0*R0 + z*z) -2.5*Tdisk,0.0) );

                  real drho = (Vc(RHO,k,j,i)-rhoTarget) / tauVel;
                  real dmx1 = Vc(RHO,k,j,i)*Vc(VX1,k,j,i) / tauVel + Vc(VX1,k,j,i) * drho;
                  real dmx2 = Vc(RHO,k,j,i)*Vc(VX2,k,j,i) / tauVel + Vc(VX2,k,j,i) * drho;
                  real dmx3 = Vc(RHO,k,j,i)*(Vc(VX3,k,j,i)-vx3Target) / tauVel + Vc(VX3,k,j,i) * drho;
                  real deng = Vc(VX1,k,j,i)*dmx1 + Vc(VX2,k,j,i)*dmx2 + Vc(VX3,k,j,i)*dmx3;

                  Uc(RHO,k,j,i) += -drho*dt;
                  Uc(MX1,k,j,i) += -dmx1*dt;
                  Uc(MX2,k,j,i) += -dmx2*dt;
                  Uc(MX3,k,j,i) += -dmx3*dt;
                  Uc(ENG,k,j,i) += -deng*dt;
                }*/

});


}




void InternalBoundary(Fluid<DefaultPhysics> *hydro, const real t) {
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray4D<real> dustVc;
  auto *data = hydro->data;
  bool haveDust = data->haveDust;
  if(haveDust) {
   dustVc= data->dust[0]->Vc;
  }
  IdefixArray4D<real> Vs = data->hydro->Vs;
  IdefixArray1D<real> x1=data->x[IDIR];
  IdefixArray1D<real> x2=data->x[JDIR];

  real vAmax = 2.0;
  real densityFloor0 = densityFloorGlob;
  real Rin = 1.0;
  real epsilon=epsilonGlob;

  idefix_for("InternalBoundary",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
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
      real densityFloor = computeDensityFloor(R,z,densityFloor0,Rin,epsilon);
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
        real Omega=1.0;
        real Rin = 1.0;
        real csdisk = epsilonGlob/sqrt(Rin);
        real cscorona = epsilonTopGlob/sqrt(Rin);
        real densityFloor0 = densityFloorGlob;
        real epsilon=epsilonGlob;

        hydro->boundary->BoundaryFor("UserDefX1",dir,side,
            KOKKOS_LAMBDA (int k, int j, int i) {
                real R=x1(i)*sin(x2(j));
                real z=x1(i)*cos(x2(j));
                /*
                Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost);
                Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost);*/

                //Vc(RHO,k,j,i) = 1.0/(Rin*sqrt(Rin))  * exp(1.0/ (csdisk*csdisk) * (1.0/sqrt(Rin*Rin+z*z)-1.0/Rin));
                //real densityFloor = computeDensityFloor(R,z,densityFloor0,Rin,epsilon);
                //if(Vc(RHO,k,j,i) < densityFloor) Vc(RHO,k,j,i) = densityFloor;

                Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost);

                Vc(PRS,k,j,i) = Vc(RHO,k,j,i)*csdisk*csdisk;

                if(Vc(VX1,k,j,ighost)>=ZERO_F) Vc(VX1,k,j,i) = -Vc(VX1,k,j,2*ighost-i-1);
                       else Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost);
                Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost);
                //real Rmin = FMAX(0.3,R);

                //Vc(VX3,k,j,i) = 1.0/sqrt(Rmin) * sqrt( Rmin / sqrt(Rmin*Rmin + z*z));
                Vc(VX3,k,j,i) = Omega*R;
                #if DIMENSIONS < 3
                Vc(BX3,k,j,i) = - Vc(BX3,k,j,2*ighost-i-1);
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
        real Rin = 1.0;
        real csdisk = epsilonGlob/sqrt(Rin);
        real cscorona = epsilonTopGlob/sqrt(Rin);

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


void EmfBoundary(Fluid<DefaultPhysics> *hydro, const real t) {
    IdefixArray3D<real> Ex1 = hydro->emf->ex;
    IdefixArray3D<real> Ex2 = hydro->emf->ey;
    IdefixArray3D<real> Ex3 = hydro->emf->ez;
    auto *data = hydro->data;

    if(data->lbound[IDIR] == userdef) {

        int ighost = data->nghost[IDIR];

        // Do not permit poloidal field to enter the seed
        idefix_for("EMFBoundary",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,ighost+1,
                    KOKKOS_LAMBDA (int k, int j, int i) {
            Ex3(k,j,i) = ZERO_F;
            #if DIMENSIONS == 3
              Ex2(k,j,i) = ZERO_F;
            #endif
        });
    }
}

void FluxBoundary(DataBlock & data, int dir, BoundarySide side, const real t) {
    IdefixArray4D<real> Flux = data.hydro->FluxRiemann;
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
        int c = 1.0/sin(th(j));
        if(c>5) c = 5;
        coarseningLevel(j,i) = c;

      });
}


void ComputeUserVars(DataBlock & data, UserDefVariablesContainer &variables) {

  // custom scratch array
  IdefixArray3D<real> scrh("Scratch", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);

  // Mirror data on Host
  DataBlockHost d(data);

  // Sync it
  d.SyncFromDevice();

  // Make references to the user-defined arrays (variables is a container of IdefixHostArray3D)
  // Note that the labels should match the variable names in the input file
  //IdefixHostArray3D<real> VdX1  = variables["VdX1"];
  //IdefixHostArray3D<real> VdX2  = variables["VdX2"];
  //IdefixHostArray3D<real> VdX3  = variables["VdX3"];
  IdefixHostArray3D<real> InvDt  = variables["InvDt"];

  IdefixHostArray1D<real> x1=d.x[IDIR];
  IdefixHostArray1D<real> x2=d.x[JDIR];
  IdefixHostArray4D<real> Vc=d.Vc;
  IdefixHostArray4D<real> J=d.J;
  IdefixArray3D<real>::host_mirror_type scrhHost = Kokkos::create_mirror_view(scrh);
  Kokkos::deep_copy(scrhHost,scrh);
  IdefixArray3D<real>::host_mirror_type xHHost;

  for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
    for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
      for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
        real z=x1(i)*cos(x2(j));
        real R=FMAX(FABS(x1(i)*sin(x2(j))),ONE_F);
        real H=R*epsilonGlob;
        real Omega=pow(R,-1.5);
        InvDt(k,j,i) = d.InvDt(k,j,i);

        // Compute ion drift speed JxB
        // Compute J at cell center
        /*
        #if DIMENSIONS < 3
        real Jx1 = AVERAGE_4D_Y(J,IDIR,k,j+1,i);
        real Jx2 = AVERAGE_4D_X(J,JDIR,k,j,i+1);
        real Jx3 = AVERAGE_4D_XY(J,KDIR,k,j+1,i+1);
        #else
        real Jx1 = AVERAGE_4D_YZ(J,IDIR,k+1,j+1,i);
        real Jx2 = AVERAGE_4D_XZ(J,JDIR,k+1,j,i+1);
        real Jx3 = AVERAGE_4D_XY(J,KDIR,k,j+1,i+1);
        #endif
        real Bx1 = Vc(BX1,k,j,i);
        real Bx2 = Vc(BX2,k,j,i);
        real Bx3 = Vc(BX3,k,j,i);

        VdX1(k,j,i) = scrhHost(k,j,i) * (Jx2 * Bx3 - Jx3 * Bx2);
        VdX2(k,j,i) = scrhHost(k,j,i) * (Jx3 * Bx1 - Jx1 * Bx3);
        VdX3(k,j,i) = scrhHost(k,j,i) * (Jx1 * Bx2 - Jx2 * Bx1);
*/
      }
    }
  }

  #ifdef EVOLVE_VECTOR_POTENTIAL
  IdefixHostArray3D<real> psi;
  try {
    // Try to get the array
    psi  = variables.at("psi");
  } catch(std::exception &e) {
    // Array is not defined
    return;
  }
  for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
    for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
      for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
        psi(k,j,i) = 0.25*(
              d.Ve(AX3e,k,j,i)*d.xl[IDIR](i)*sin(d.xl[JDIR](j)) +
              d.Ve(AX3e,k,j,i+1)*d.xl[IDIR](i+1)*sin(d.xl[JDIR](j)) +
              d.Ve(AX3e,k,j+1,i)*d.xl[IDIR](i)*sin(d.xl[JDIR](j+1)) +
              d.Ve(AX3e,k,j+1,i+1)*d.xl[IDIR](i+1)*sin(d.xl[JDIR](j+1))
            );
  }}}

  IdefixHostArray3D<real> EphiIdeal;
  try {
    // Try to get the array
    EphiIdeal  = variables.at("EphiIdeal");
  } catch(std::exception &e) {
    // Array is not defined
    return;
  }
  // Talk to the emf object to recompute required EMFs
  IdefixArray3D<real> ex1,ex2,ex3;

  ex1 = emf->ex;
  ex2 = emf->ey;
  ex3 = emf->ez;

  // Compute Ideal EMF
  data.SetBoundaries();
  emf->CalcCornerEMF(data.t);
  Kokkos::deep_copy(EphiIdeal,ex3);

  idefix_for("resetEMF", 0, data.np_tot[KDIR], 0, data.np_tot[JDIR], 0, data.np_tot[IDIR],
      KOKKOS_LAMBDA(int k, int j,int i) {
        #if DIMENSIONS == 3
        ex1(k,j,i) = 0.0;
        ex2(k,j,i) = 0.0;
        #endif
        ex3(k,j,i) = 0.0;
      });


  // Done
  #endif
}

void analysisFunction(DataBlock& data) {
  analysis->PerformAnalysis(data);
}


void ComputeRho(DataBlock &data) {
  real Rin=1.0;
  real alpha = alphaGlob; // Density power law

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
    real Tdisk = epsilonGlob*epsilonGlob;
    real Tcorona = epsilonTopGlob*epsilonTopGlob;
    //if(x1(i) < 1.5) cscorona = csdisk;
    real Teff=0.5*(Tdisk+Tcorona)+0.5*(Tcorona-Tdisk)*tanh((Zh-HidealGlob)/trSmoothingGlob);
    // BS17 eq. 14
    logf(j) = logf(j-1)+dth*(cos(th)/sin(th)*1/Teff-(alpha+1));
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

      real Tdisk = epsilonGlob*epsilonGlob;
      real Tcorona = epsilonTopGlob*epsilonTopGlob;
      //if(x1(i) < 1.5) cscorona = csdisk;
      real Teff=0.5*(Tdisk+Tcorona)+0.5*(Tcorona-Tdisk)*tanh((Zh-HidealGlob)/trSmoothingGlob);

      // BS17 F(theta)=f(theta)*g(theta)
      // here logF is log(F) and g is Teff
      // We normalise f so that it is=1 in the disc midaplne
      real f = exp(logf(jglob)) *epsilonGlob*epsilonGlob / Teff;

      rhoH(j,i) = pow(Ri,-alpha) * f;
      real densityFloor = computeDensityFloor(R,z,densityFloorGlob,Rin,epsilonGlob);
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
  if(data.haveGridCoarsening) {
    data.EnrollGridCoarseningLevels(&CoarsenFunction);
  }
  data.hydro->EnrollUserSourceTerm(&MySourceTerm);
  data.hydro->EnrollInternalBoundary(&InternalBoundary);
  //data.hydro->EnrollEmfBoundary(&EmfBoundary);
  //data.hydro->EnrollFluxBoundary(&FluxBoundary);

  gammaGlob=data.hydro->eos->GetGamma();
  alphaGlob = input.GetOrSet<real>("Setup","alpha",0,1.5);  // density power law
  alphaBetaGlob = input.GetOrSet<real>("Setup","alphaBeta",0,0.0); // beta power law
  epsilonGlob = input.Get<real>("Setup","epsilon",0);
  epsilonTopGlob = input.Get<real>("Setup","epsilonTop",0);
  betaGlob = input.Get<real>("Setup","beta",0);
  HidealGlob = input.Get<real>("Setup","Hideal",0);

  densityFloorGlob = input.Get<real>("Setup","densityFloor",0);
  trSmoothingGlob = input.Get<real>("Setup","transitionSmoothing",0);

  output.EnrollUserDefVariables(&ComputeUserVars);
  // assume disc surface at 6.5 h
  analysis = new Analysis(grid, data,std::string("profile.dat"), HidealGlob*epsilonGlob);
  output.EnrollAnalysis(&analysisFunction);


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

    // Make vector potential
    IdefixHostArray4D<real> A = IdefixHostArray4D<real>("Setup_VectorPotential", 3, data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);

    // Get the equilibrium density profile
    IdefixHostArray2D<real> rhoH("RhoHost",d.np_tot[JDIR],d.np_tot[IDIR]);
    Kokkos::deep_copy(rhoH,*rhoInit);

    real Rin=1.0;
    real alpha = alphaGlob; // MINUS Density power law
    real gammaB =0.5*(-alpha-1-alphaBetaGlob); // MAgnetic field power law

    real B0 = epsilonGlob*sqrt(2.0/betaGlob);

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


                real Tdisk = epsilonGlob*epsilonGlob;
                real Tcorona = epsilonTopGlob*epsilonTopGlob;
                //if(x1(i) < 1.5) cscorona = csdisk;
                real Teff=0.5*(Tdisk+Tcorona)+0.5*(Tcorona-Tdisk)*tanh((Zh-HidealGlob)/trSmoothingGlob);

                d.Vc(RHO,k,j,i) = rhoH(j,i);
                d.Vc(PRS,k,j,i) = Teff/Ri * d.Vc(RHO,k,j,i);
                // BS17 eq. 15
                d.Vc(VX3,k,j,i) = sqrt((1-(alpha+1)*Teff)/Ri);

                // Enforce solid body rotation above the "seed"
                if(R<Rin) d.Vc(VX3,k,j,i) = sqrt((1-(alpha+1)*Teff))*R/Rin;

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
                                      1/(gammaB+2)*(pow(R,gammaB+1) - pow(Rin,gammaB+2)/R)
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
