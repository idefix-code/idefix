// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file test/Radiation/novae_like_3D_no_damping_rad/setup.cpp
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
#include "dumpImage.hpp"

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
real ToutGlob;
real TceilingGlob;
real TfloorGlob;
real kappapGlob;
real kapparGlob;
std::string kappatypeGlob;
real kramersTindexGlob;
real kramersrhoindexGlob;

LookupTable<2> *kappaptabGlob;
LookupTable<2> *kappartabGlob;

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


void InternalBoundary(Fluid<DefaultPhysics> *hydro, const real t) {
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray4D<real> Uc = hydro->Uc;
  auto *data = hydro->data;
  IdefixArray4D<real> Vs = data->hydro->Vs;
  IdefixArray1D<real> x1=data->x[IDIR];
  IdefixArray1D<real> x2=data->x[JDIR];

  real vAmax = 2.0;
  real densityFloor0 = densityFloorGlob;
  real Rin = 1.0;
  real epsilon=epsilonGlob;

  auto units = idfx::units;
  real mu = muGlob;
  real gamma = gammaGlob;
  real Tceiling = TceilingGlob;
  real Tfloor = TfloorGlob;

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
      real T = Vc(PRS,k,j,i)/Vc(RHO,k,j,i)*units.GetKelvin()*mu;
      real Eint_old = Vc(PRS,k,j,i)/(gamma-1.0);
      if (T<Tfloor){
        Vc(PRS,k,j,i)=Tfloor*Vc(RHO,k,j,i)/(units.GetKelvin()*mu);
        Uc(ENG,k,j,i) -= Eint_old;
        Uc(ENG,k,j,i) += Vc(PRS,k,j,i)/(gamma-1.0);
      } else if (T>Tceiling){
        Vc(PRS,k,j,i)=Tceiling*Vc(RHO,k,j,i)/(units.GetKelvin()*mu);
        Uc(ENG,k,j,i) -= Eint_old;
        Uc(ENG,k,j,i) += Vc(PRS,k,j,i)/(gamma-1.0);
      }

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

        hydro->boundary->BoundaryFor("UserDefX1",dir,side,
            KOKKOS_LAMBDA (int k, int j, int i) {
                real R=x1(i)*sin(x2(j));
                real z=x1(i)*cos(x2(j));

                Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost);

                Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost);

                if(Vc(VX1,k,j,ighost)>=ZERO_F) Vc(VX1,k,j,i) = 0.0;
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
                Vc(VX3,k,j,i) = Vc(VX3,k,j,ighost);
                #if DIMENSIONS < 3
                Vc(BX3,k,j,i) = - Vc(BX3,k,j,2*ighost-i+1);
                #endif

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

    idefix_for("EMFBoundary",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
            Ex1(k,j,i) = ZERO_F;
            Ex2(k,j,i) = ZERO_F;
            Ex3(k,j,i) = ZERO_F;
        });
}

void FluxBoundary(DataBlock & data, int dir, BoundarySide side, const real t) {
    IdefixArray4D<real> Flux = data.hydro->FluxRiemann[dir];


    idefix_for("FluxInternal",
                0, data.np_tot[KDIR],
                0, data.np_tot[JDIR],
                0, data.np_tot[IDIR],
       KOKKOS_LAMBDA (int k, int j, int i) {
         Flux(RHO, k, j, i) = 0.0;
         Flux(MX1, k, j, i) = 0.0;
         Flux(MX2, k, j, i) = 0.0;
         Flux(MX3, k, j, i) = 0.0;
         Flux(ENG, k, j, i) = 0.0;
     });
}

void FluxBoundaryRad(DataBlock & data, int dir, BoundarySide side, const real t) {
   IdefixArray4D<real> Flux = data.radiation[0]->FluxRiemann[dir];

  if(dir==IDIR) {
    int ighost,nxi,iend,ibeg;
    if(side == left) {
      ighost = data.nghost[IDIR];
      ibeg = 0;
      iend = data.beg[IDIR];
      idefix_for("UserDefBoundaryRad",
        0, data.np_tot[KDIR],
        0, data.np_tot[JDIR],
        ibeg, iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          if (Flux(ER,k,j,i) >= ZERO_F) Flux(ER,k,j,i) = ZERO_F;
          if (Flux(FR1,k,j,i) >= ZERO_F) Flux(FR1,k,j,i) = ZERO_F;
          if (Flux(FR2,k,j,i) >= ZERO_F) Flux(FR2,k,j,i) = ZERO_F;
        });
    } else if (side==right){
      ighost = data.nghost[IDIR];
      nxi = data.np_int[IDIR];
      ibeg = data.end[IDIR];
      iend =data.np_tot[IDIR];
      idefix_for("UserDefBoundaryRad",
        0, data.np_tot[KDIR],
        0, data.np_tot[JDIR],
        ibeg, iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          if (Flux(ER,k,j,i) <= ZERO_F) Flux(ER,k,j,i) = ZERO_F;
          if (Flux(FR1,k,j,i) <= ZERO_F) Flux(FR1,k,j,i) = ZERO_F;
          if (Flux(FR2,k,j,i) <= ZERO_F) Flux(FR2,k,j,i) = ZERO_F;
        });
    }
   }

   if(dir==JDIR) {
    int jghost,nxj,jend,jbeg;
    if(side == left) {
      jghost = data.nghost[JDIR];
      jbeg = 0;
      jend = data.beg[JDIR];
      idefix_for("UserDefBoundaryRad",
        0, data.np_tot[KDIR],
        jbeg, jend,
        0, data.np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
          if (Flux(ER,k,j,i) >= ZERO_F) Flux(ER,k,j,i) = ZERO_F;
          if (Flux(FR1,k,j,i) >= ZERO_F) Flux(FR1,k,j,i) = ZERO_F;
          if (Flux(FR2,k,j,i) >= ZERO_F) Flux(FR2,k,j,i) = ZERO_F;
        });
    } else if (side==right){
      jghost = data.nghost[JDIR];
      nxj = data.np_int[JDIR];
      jbeg = data.end[JDIR];
      jend =data.np_tot[JDIR];
      idefix_for("UserDefBoundaryRad",
        0, data.np_tot[KDIR],
        jbeg, jend,
        0, data.np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
          if (Flux(ER,k,j,i) <= ZERO_F) Flux(ER,k,j,i) = ZERO_F;
          if (Flux(FR1,k,j,i) <= ZERO_F) Flux(FR1,k,j,i) = ZERO_F;
          if (Flux(FR2,k,j,i) <= ZERO_F) Flux(FR2,k,j,i) = ZERO_F;
        });
    }
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


void InternalBoundaryRad(Fluid<RadiationPhysics> *radiation, const real t) {
  IdefixArray4D<real> Vc = radiation->Vc;
  auto *data = radiation->data;
  auto units = idfx::units;

  real Tceiling = TceilingGlob;

  idefix_for("InternalBoundaryRad",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
        if (Vc(ER,k,j,i) > units.ar*std::pow(Tceiling,4)/units.GetEnergy()){
          Vc(ER,k,j,i) = units.ar*std::pow(Tceiling,4)/units.GetEnergy();
          real Fnorm = EXPAND(Vc(FR1,k,j,i)*Vc(FR1,k,j,i),+Vc(FR2,k,j,i)*Vc(FR2,k,j,i),+Vc(FR3,k,j,i)*Vc(FR3,k,j,i));
          if (Fnorm > Vc(ER,k,j,i)){
            EXPAND(Vc(FR1,k,j,i) *= Vc(ER,k,j,i)/Fnorm;,
                   Vc(FR2,k,j,i) *= Vc(ER,k,j,i)/Fnorm;,
                   Vc(FR3,k,j,i) *= Vc(ER,k,j,i)/Fnorm;)
          }
        }
        //std::printf("Er=%e in InternalBoundary\n",Vc(ER,k,j,i));
    });
}

void UserdefBoundaryRad(Fluid<RadiationPhysics> *radiation, int dir, BoundarySide side, real t) {
  IdefixArray4D<real> Vc = radiation->Vc;
  auto *data = radiation->data;
  auto units=idfx::units;
  IdefixArray4D<real> FluxRad = radiation->FluxRiemann[dir];


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
          Vc(ER,k,j,i) = Vc(ER,k,j,ighost);
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
          Vc(ER,k,j,i) = Vc(ER,k,j,ighost+nxi-1);
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

  if(dir==JDIR) {
    int jghost,nxj,jend,jbeg;
    if(side == left) {
      jghost = data->nghost[JDIR];
      jbeg = 0;
      jend = data->beg[JDIR];
      idefix_for("UserDefBoundaryRad",
        0, data->np_tot[KDIR],
        jbeg, jend,
        0, data->np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vc(ER,k,j,i) = Vc(ER,k,jghost,i);
          if (Vc(FR2,k,jghost,i) >=ZERO_F){
            Vc(FR2,k,j,i) = ZERO_F;
          } else {
            Vc(FR2,k,j,i) = Vc(FR2,k,jghost,i);
          }
          Vc(FR1,k,j,i) = Vc(FR1,k,jghost,i);
          Vc(FR3,k,j,i) = Vc(FR3,k,jghost,i);
        });
    } else if (side==right){
      jghost = data->nghost[JDIR];
      nxj = data->np_int[JDIR];
      jbeg = data->end[JDIR];
      jend =data->np_tot[JDIR];
      idefix_for("UserDefBoundaryRad",
        0, data->np_tot[KDIR],
        jbeg, jend,
        0, data->np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vc(ER,k,j,i) = Vc(ER,k,jghost+nxj-1,i);
          if (Vc(FR2,k,jghost+nxj-1,i) <=ZERO_F){
            Vc(FR2,k,j,i) = ZERO_F;
            //FluxRad(ER,k,j,i) = ZERO_F;
          } else {
            Vc(FR2,k,j,i) = Vc(FR2,k,jghost+nxj-1,i);
          }
          Vc(FR1,k,j,i) = Vc(FR1,k,jghost+nxj-1,i);
          Vc(FR3,k,j,i) = Vc(FR3,k,jghost+nxj-1,i);
        });
    }
  }

}

void ComputeUserVars(DataBlock & data, UserDefVariablesContainer &variables) {

  // custom scratch array
  IdefixArray3D<real> scrh1("Scratch1", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);
  IdefixArray3D<real> scrh2("Scratch2", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);
  IdefixArray3D<real> scrh3("Scratch3", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);
  IdefixArray3D<real> scrh4("Scratch4", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);
  IdefixArray3D<real> scrh5("Scratch5", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);

  IdefixArray4D<real> FluxRiemannIDIR = data.hydro->FluxRiemann[IDIR];
  IdefixArray4D<real> FluxRiemannJDIR = data.hydro->FluxRiemann[JDIR];
  IdefixArray4D<real> RadFluxRiemannIDIR = data.radiation[0]->FluxRiemann[IDIR];
  IdefixArray4D<real> RadFluxRiemannJDIR = data.radiation[0]->FluxRiemann[JDIR];

  idefix_for("UserVar",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
   KOKKOS_LAMBDA (int k, int j, int i) {
      scrh1(k,j,i) = FluxRiemannIDIR(RHO,k,j,i);
      scrh2(k,j,i) = FluxRiemannJDIR(RHO,k,j,i);
      scrh3(k,j,i) = FluxRiemannIDIR(VX1,k,j,i);
      scrh4(k,j,i) = FluxRiemannJDIR(VX1,k,j,i);
      scrh5(k,j,i) = FluxRiemannIDIR(VX2,k,j,i);
   });
  Kokkos::deep_copy(variables["Fluxrhor"], scrh1);
  Kokkos::deep_copy(variables["Fluxrhot"], scrh2);
  Kokkos::deep_copy(variables["Fluxmrr"], scrh3);
  Kokkos::deep_copy(variables["Fluxmrt"], scrh4);
  Kokkos::deep_copy(variables["Fluxmtr"], scrh5);

  idefix_for("UserVar",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
   KOKKOS_LAMBDA (int k, int j, int i) {
     scrh1(k,j,i) = FluxRiemannJDIR(VX2,k,j,i);
     scrh2(k,j,i) = FluxRiemannIDIR(ENG,k,j,i);
     scrh3(k,j,i) = FluxRiemannJDIR(ENG,k,j,i);
     scrh4(k,j,i) = RadFluxRiemannIDIR(ER,k,j,i);
     scrh5(k,j,i) = RadFluxRiemannJDIR(ER,k,j,i);
   });
  Kokkos::deep_copy(variables["Fluxmtt"], scrh1);
  Kokkos::deep_copy(variables["FluxEngr"], scrh2);
  Kokkos::deep_copy(variables["FluxEngt"], scrh3);
  Kokkos::deep_copy(variables["FluxRadErr"], scrh4);
  Kokkos::deep_copy(variables["FluxRadErt"], scrh5);

  idefix_for("UserVar",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
   KOKKOS_LAMBDA (int k, int j, int i) {
     scrh1(k,j,i) = RadFluxRiemannIDIR(FR1,k,j,i);
     scrh2(k,j,i) = RadFluxRiemannJDIR(FR1,k,j,i);
     scrh3(k,j,i) = RadFluxRiemannIDIR(FR2,k,j,i);
     scrh4(k,j,i) = RadFluxRiemannJDIR(FR2,k,j,i);
   });
  Kokkos::deep_copy(variables["FluxRadFrr"], scrh1);
  Kokkos::deep_copy(variables["FluxRadFrt"], scrh2);
  Kokkos::deep_copy(variables["FluxRadFtr"], scrh3);
  Kokkos::deep_copy(variables["FluxRadFtt"], scrh4);

  // Mirror data on Host
  DataBlockHost d(data);

  // Sync it
  d.SyncFromDevice();

  auto units = idfx::units;

  // Make references to the user-defined arrays (variables is a container of IdefixHostArray3D)
  // Note that the labels should match the variable names in the input file

  for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
    for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
      for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
        variables["rhovrvr"](k,j,i) = d.Vc(RHO,k,j,i)*d.Vc(VX1,k,j,i)*d.Vc(VX1,k,j,i);
        variables["rhovrvt"](k,j,i) = d.Vc(RHO,k,j,i)*d.Vc(VX1,k,j,i)*d.Vc(VX2,k,j,i);
        variables["rhovrvp"](k,j,i) = d.Vc(RHO,k,j,i)*d.Vc(VX1,k,j,i)*d.Vc(VX3,k,j,i);
        variables["rhovtvt"](k,j,i) = d.Vc(RHO,k,j,i)*d.Vc(VX2,k,j,i)*d.Vc(VX2,k,j,i);
        variables["rhovtvp"](k,j,i) = d.Vc(RHO,k,j,i)*d.Vc(VX2,k,j,i)*d.Vc(VX3,k,j,i);
        variables["rhovpvp"](k,j,i) = d.Vc(RHO,k,j,i)*d.Vc(VX3,k,j,i)*d.Vc(VX3,k,j,i);
        variables["BrBr"](k,j,i) = d.Vc(BX1,k,j,i)*d.Vc(BX1,k,j,i);
        variables["BrBt"](k,j,i) = d.Vc(BX1,k,j,i)*d.Vc(BX2,k,j,i);
        variables["BrBp"](k,j,i) = d.Vc(BX1,k,j,i)*d.Vc(BX3,k,j,i);
        variables["BtBt"](k,j,i) = d.Vc(BX2,k,j,i)*d.Vc(BX2,k,j,i);
        variables["BtBp"](k,j,i) = d.Vc(BX2,k,j,i)*d.Vc(BX3,k,j,i);
        variables["BpBp"](k,j,i) = d.Vc(BX3,k,j,i)*d.Vc(BX3,k,j,i);
        real rhov2 = d.Vc(RHO,k,j,i)*(d.Vc(VX1,k,j,i)*d.Vc(VX1,k,j,i)+d.Vc(VX2,k,j,i)*d.Vc(VX2,k,j,i)+d.Vc(VX3,k,j,i)*d.Vc(VX3,k,j,i));
        variables["rhov2vr"](k,j,i) = rhov2*d.Vc(VX1,k,j,i);
        variables["rhov2vt"](k,j,i) = rhov2*d.Vc(VX2,k,j,i);
        real B2 = d.Vc(BX1,k,j,i)*d.Vc(BX1,k,j,i)+d.Vc(BX2,k,j,i)*d.Vc(BX2,k,j,i)+d.Vc(BX3,k,j,i)*d.Vc(BX3,k,j,i);
        variables["B2vr"](k,j,i) = B2*d.Vc(VX1,k,j,i);
        variables["B2vt"](k,j,i) = B2*d.Vc(VX2,k,j,i);
        real BV = d.Vc(BX1,k,j,i)*d.Vc(VX1,k,j,i)+d.Vc(BX2,k,j,i)*d.Vc(VX2,k,j,i)+d.Vc(BX3,k,j,i)*d.Vc(VX3,k,j,i);
        variables["BVBr"](k,j,i) = BV*d.Vc(BX1,k,j,i);
        variables["BVBt"](k,j,i) = BV*d.Vc(BX2,k,j,i);
        variables["Emfr"](k,j,i) = d.Ex1(k,j,i);
        variables["Emft"](k,j,i) = d.Ex2(k,j,i);
        variables["Emfp"](k,j,i) = d.Ex3(k,j,i);
        real T = d.Vc(PRS,k,j,i)/d.Vc(RHO,k,j,i)*units.GetKelvin()*muGlob;
        real logrho = std::log10(d.Vc(RHO,k,j,i)*units.GetDensity());
        real x[2];
        x[0] = FMIN(-4.05,FMAX(-14.,logrho));
        x[1] = FMIN(FMAX(std::log10(T),2.5),5.98);
        if (kappatypeGlob.compare("usertable") == 0){
          variables["kappap"](k,j,i) = std::pow(10.,kappaptabGlob->GetHost(x));
          variables["kappar"](k,j,i) = std::pow(10.,kappartabGlob->GetHost(x));
        } else if (kappatypeGlob.compare("constant") == 0) {
          variables["kappap"](k,j,i) = kappapGlob;
          variables["kappar"](k,j,i) = kapparGlob;
        } else if (kappatypeGlob.compare("kramers") == 0) {
          variables["kappap"](k,j,i) = kappapGlob*d.Vc(RHO,k,j,i)*units.GetDensity()/kramersrhoindexGlob*std::pow(T/kramersTindexGlob,-3.5);
          variables["kappar"](k,j,i) = kapparGlob*d.Vc(RHO,k,j,i)*units.GetDensity()/kramersrhoindexGlob*std::pow(T/kramersTindexGlob,-3.5);
        }
        variables["rhokappapEr"](k,j,i) = d.Vc(RHO,k,j,i)*units.GetDensity()*variables["kappap"](k,j,i)*d.RadVc[0](ER,k,j,i)*units.GetEnergy();
        variables["rhokapparFr"](k,j,i) = d.Vc(RHO,k,j,i)*units.GetDensity()*variables["kappar"](k,j,i)*d.RadVc[0](FR1,k,j,i)*units.GetEnergy();
        variables["rhokapparFt"](k,j,i) = d.Vc(RHO,k,j,i)*units.GetDensity()*variables["kappar"](k,j,i)*d.RadVc[0](FR2,k,j,i)*units.GetEnergy();
        variables["rhokappaparT4"](k,j,i) = d.Vc(RHO,k,j,i)*units.GetDensity()*variables["kappap"](k,j,i)*units.ar*std::pow(T,4.);
        variables["T"](k,j,i) = T;
      }
    }
  }
}

void analysisFunction(DataBlock& data) {
  analysis->PerformAnalysis(data);
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
  data.hydro->EnrollInternalBoundary(&InternalBoundary);
  //data.hydro->EnrollEmfBoundary(&EmfBoundary);
  //data.hydro->EnrollFluxBoundary(&FluxBoundary);

  if(data.haveRadiation) {
    int nFrequencies = data.radiation.size();
    data.radiation[0]->EnrollUserDefBoundary(&UserdefBoundaryRad);
    data.radiation[0]->EnrollInternalBoundary(&InternalBoundaryRad);
    //data.radiation[0]->EnrollFluxBoundary(&FluxBoundaryRad);
  }

  gammaGlob=data.hydro->eos->GetGamma();
  muGlob = input.Get<real>("Hydro","mu",0);
  alphaGlob = input.GetOrSet<real>("Setup","alpha",0,1.5);  // density power law
  alphaBetaGlob = input.GetOrSet<real>("Setup","alphaBeta",0,0.0); // beta power law
  epsilonGlob = input.Get<real>("Setup","epsilon",0);
  epsilonTopGlob = input.Get<real>("Setup","epsilonTop",0);
  betaGlob = input.Get<real>("Setup","beta",0);
  HidealGlob = input.Get<real>("Setup","Hideal",0);
  TceilingGlob = input.Get<real>("Setup","Tceiling",0);
  TfloorGlob = input.Get<real>("Setup","Tfloor",0);

  densityFloorGlob = input.Get<real>("Setup","densityFloor",0);
  trSmoothingGlob = input.Get<real>("Setup","transitionSmoothing",0);

  output.EnrollUserDefVariables(&ComputeUserVars);
  // assume disc surface at 6.5 h
  analysis = new Analysis(grid, data,std::string("profile.dat"), HidealGlob*epsilonGlob);
  output.EnrollAnalysis(&analysisFunction);

  kappatypeGlob = input.Get<std::string>("Rad","kappa",0);
  if (kappatypeGlob.compare("usertable") == 0){
    std::string kappapfile = input.Get<std::string>("Rad","kappa",2);
    std::string kapparfile = input.Get<std::string>("Rad","kappa",3);
    kappaptabGlob = new LookupTable<2>(kappapfile,',');
    kappartabGlob = new LookupTable<2>(kapparfile,',');
   } else if (kappatypeGlob.compare("constant") == 0) {
    kappapGlob = input.Get<real>("Rad","kappa",1);
    kapparGlob = input.Get<real>("Rad","kappa",2);
   } else if (kappatypeGlob.compare("kramers") == 0) {
    kramersTindexGlob = input.Get<real>("Rad","kappa",4);
    kramersrhoindexGlob = input.Get<real>("Rad","kappa",3);
    kappapGlob = input.Get<real>("Rad","kappa",1);
    kapparGlob = input.Get<real>("Rad","kappa",2);
   }
}


// Flow initialisation, read directly from the DumpImage
void Setup::InitFlow(DataBlock &data) {

  // Create a host copy
  DataBlockHost d(data);

  DumpImage image("dump.0241.dmp", &data);

  for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
    for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
      for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {

        // Note that the restart dump array only contains the full (global) active domain
        // (i.e. it excludes the boundaries, but it is not decomposed accross MPI procs)
        int iglob=i-2*d.beg[IDIR]+d.gbeg[IDIR];
        int jglob=j-2*d.beg[JDIR]+d.gbeg[JDIR];
        int kglob=k-2*d.beg[KDIR]+d.gbeg[KDIR];

        d.Vc(RHO,k,j,i) = image.arrays["Vc-RHO"](kglob,jglob,iglob);
        d.Vc(PRS,k,j,i) = image.arrays["Vc-PRS"](kglob,jglob,iglob);
        real T = d.Vc(PRS,k,j,i)/d.Vc(RHO,k,j,i)*idfx::units.GetKelvin()*muGlob;
        d.RadVc[0](ER,k,j,i) = idfx::units.ar*std::pow(T,4)/idfx::units.GetEnergy();
        d.RadVc[0](FR1,k,j,i) = d.RadVc[0](ER,k,j,i);
        d.RadVc[0](FR2,k,j,i) = 0.;
        d.RadVc[0](FR3,k,j,i) = 0.;
        d.Vc(VX1,k,j,i) = image.arrays["Vc-VX1"](kglob,jglob,iglob);
        d.Vc(VX2,k,j,i) = image.arrays["Vc-VX2"](kglob,jglob,iglob);
        d.Vc(VX3,k,j,i) = image.arrays["Vc-VX3"](kglob,jglob,iglob);
}}}

  // For magnetic variable, we should fill the entire active domain, hence an additional
  // point in the field direction
  for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
    for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR]+IOFFSET ; i++) {
          int iglob=i-2*d.beg[IDIR]+d.gbeg[IDIR];
          int jglob=j-2*d.beg[JDIR]+d.gbeg[JDIR];
          int kglob=k-2*d.beg[KDIR]+d.gbeg[KDIR];
          d.Vs(BX1s,k,j,i) = image.arrays["Vs-BX1s"](kglob,jglob,iglob);
          d.Vs(BX2s,k,j,i) = image.arrays["Vs-BX2s"](kglob,jglob,iglob);
          d.Vs(BX3s,k,j,i) = image.arrays["Vs-BX3s"](kglob,jglob,iglob);
  }}}


  // Send our datablock to the device
  d.SyncToDevice();
}




Setup::~Setup() {
  delete analysis;
}
