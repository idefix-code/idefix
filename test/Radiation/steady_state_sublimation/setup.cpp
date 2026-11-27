// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file test/Radiation/steady_state_sublimation/setup.cpp
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
#include "column.hpp"
#include "lookupTable.hpp"

real R0Glob;
real rho0Glob;
real rhominGlob;
real T0Glob;
real muGlob;
real gammaGlob;
real rsGlob;
real TsGlob;
real kappaGlob;
real kappairrGlob;
std::string kappatypeGlob;
real alphaMRIGlob;
real alphaDZGlob;
real epsilonGlob;
real TsubGlob;
real rhosGlob;
real rhosindexGlob;
real TwidthGlob;
real f0Glob;
real kappastarGlob;
real kappagasGlob;
real MGlob;
real GGlob;
real densityFloorGlob;
real TMriGlob;
real rhoindexGlob;
real alphawidthGlob;

Column *columnGlob;
Column *columnGlob2;
LookupTable<1> *kappatableGlob;


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

  columnGlob->ComputeColumn(kapparho);
  tau = columnGlob->GetColumn();

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
                //f_gtod = f0;
                kappap(k,j,i) = kappa_star*f_gtod+kappa_gas;
                kappar(k,j,i) = kappa_star*f_gtod+kappa_gas;
              });
}


void MyViscosity(DataBlock &data, const real t, IdefixArray3D<real> &eta1, IdefixArray3D<real> &eta2) {
  auto units = idfx::units;

  IdefixArray3D<real> InvDt = data.hydro->InvDt;
  IdefixArray4D<real> Vc=data.hydro->Vc;
  IdefixArray1D<real> r=data.x[IDIR];
  IdefixArray1D<real> th=data.x[JDIR];
  real alphaDZ = alphaDZGlob;
  real alphaMRI = alphaMRIGlob;
  real CG = MGlob*GGlob;
  real R0 = R0Glob;
  real T_MRI = TMriGlob;
  real mu = muGlob;
  real alphawidth = alphawidthGlob;

  idefix_for("MyViscosity",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real R = FMAX(r(i)*sin(th(j)),R0);
                real cs2 = Vc(PRS,k,j,i)/Vc(RHO,k,j,i);
                real T = cs2*units.GetKelvin()*mu;
                real Omega = std::sqrt(CG)*std::pow(R,-1.5);
                real alpha = (alphaMRI-alphaDZ)*0.5*(1.-std::tanh((T_MRI-T)/alphawidth))+alphaDZ;
                eta1(k,j,i) = alpha*cs2*Vc(RHO,k,j,i)/Omega;
                eta2(k,j,i) = 0.;
              });

}

void UserdefBoundaryNoStress(Fluid<DefaultPhysics> *hydro, int dir, BoundarySide side, real t) {
  IdefixArray4D<real> Vc = hydro->Vc;
  auto *data = hydro->data;
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  real rhomin = densityFloorGlob/idfx::units.GetDensity();

  if(dir==IDIR) {
    int ighost,nxi,iend,ibeg;
    if(side == left) {
      ighost = data->nghost[IDIR];
      ibeg = 0;
      iend = data->beg[IDIR];
      idefix_for("UserDefBoundary",
        0, data->np_tot[KDIR],
        0, data->np_tot[JDIR],
        ibeg, iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost);
          Vc(VX3,k,j,i) = x1(i)*std::sin(x2(j))*Vc(VX3,k,j,ighost)/(x1(ighost)*std::sin(x2(j)));

          Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost)/Vc(RHO,k,j,ighost)*Vc(RHO,k,j,i);
          if(Vc(VX1,k,j,ighost)>=ZERO_F){
            Vc(VX1,k,j,i)=ZERO_F;
          }else {
            Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost);
          }
          Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost);
        });
    } else if (side ==right) {
      ighost = data->nghost[IDIR];
      nxi = data->np_int[IDIR];
      ibeg = data->end[IDIR];
      iend = data->np_tot[IDIR];
      idefix_for("UserDefBoundary",
        0, data->np_tot[KDIR],
        0, data->np_tot[JDIR],
        ibeg, iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost+nxi-1);
          Vc(VX3,k,j,i) = Vc(VX3,k,j,ighost+nxi-1);

          Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost+nxi-1)/Vc(RHO,k,j,ighost+nxi-1)*Vc(RHO,k,j,i);
          if(Vc(VX1,k,j,ighost+nxi-1)<=ZERO_F){
            Vc(VX1,k,j,i)=ZERO_F;
          }else {
            Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost+nxi-1);
          }
          Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost+nxi-1);
        });
    }
  }
}

void UserdefBoundaryRad(Fluid<RadiationPhysics> *radiation, int dir, BoundarySide side, real t) {
  IdefixArray4D<real> Vc = radiation->Vc;
  auto *data = radiation->data;
  auto units=idfx::units;

  real T0 = T0Glob;
  real flux_stellar = std::pow(rsGlob/units.GetLength(),2.)*units.sigma_sb*std::pow(TsGlob,4.)/units.GetEnergy()/units.c;

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
          Vc(ER,k,j,i) = units.ar*std::pow(10.,4.)/units.GetEnergy();
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

void FluxBoundary(DataBlock & data, int dir, BoundarySide side, const real t) {
    IdefixArray4D<real> Flux = data.hydro->FluxRiemann[dir];


    idefix_for("FluxInternal",
                0, data.np_tot[KDIR],
                0, data.np_tot[JDIR],
                0, data.np_tot[IDIR],
       KOKKOS_LAMBDA (int k, int j, int i) {
         Flux(RHO, k, j, i) = 0.0;
         //Flux(MX1, k, j, i) = 0.0;
         //Flux(MX2, k, j, i) = 0.0;
         //Flux(MX3, k, j, i) = 0.0;
     });
}

// Compute user variables which will be written in vtk files
void ComputeUserVars(DataBlock & data, UserDefVariablesContainer &variables) {
  // Mirror data on Host
  DataBlockHost d(data);

  // Sync it
  d.SyncFromDevice();

  auto units = idfx::units;

  IdefixHostArray1D<real> x1=d.x[IDIR];
  IdefixHostArray1D<real> x1l=d.xl[IDIR];
  IdefixHostArray1D<real> x2=d.x[JDIR];
  IdefixHostArray1D<real> x2l=d.xl[JDIR];
  IdefixHostArray3D<real> A1=d.A[IDIR];
  IdefixHostArray3D<real> A2=d.A[JDIR];
  IdefixHostArray3D<real> dV=d.dV;

  IdefixArray3D<real> tau;
  //IdefixArray3D<real> tau2;
  IdefixArray4D<real> Vc=(&data)->hydro->Vc;
  IdefixArray3D<real> kappa=(&data)->radiation[0]->radsource->kappapArr;

  IdefixHostArray3D<real> divF  = variables["divF"];
  IdefixHostArray3D<real> A1_out  = variables["A1"];

  real rs = rsGlob;
  real Ts = TsGlob;
  real kappa_irr = kappairrGlob;

  real kirr = kappa_irr*idfx::units.GetDensity()*idfx::units.GetLength();
  real flux_pre = std::pow(rs/idfx::units.GetLength(),2.)*idfx::units.sigma_sb*std::pow(Ts,4.)/idfx::units.GetLength();


  IdefixArray3D<real> kapparho("rho",d.np_tot[KDIR],d.np_tot[JDIR],d.np_tot[IDIR]);
  idefix_for("init rho",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
    KOKKOS_LAMBDA(int k, int j, int i) {
      kapparho(k,j,i) = Vc(RHO,k,j,i)*kappa(k,j,i)*units.GetDensity()*units.GetLength();
    });

  columnGlob->ComputeColumn(kapparho);
  tau = columnGlob->GetColumn();

  //columnGlob->ComputeColumn(Vc,RHO);
  //tau = columnGlob->GetColumn();


  Kokkos::deep_copy(variables["tau"], tau);
  //Kokkos::deep_copy(variables["tau2"], tau2);
  Kokkos::deep_copy(variables["dV"], dV);

  if(kappatypeGlob=="constant") {

    for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
          A1_out(k,j,i) = A1(k,j,i);

          real Fim = std::exp(-kirr*variables["tau"](k,j,i-1))/std::pow(x1l(i),2.);
          real Fip = std::exp(-kirr*variables["tau"](k,j,i))/std::pow(x1l(i+1),2.);

          divF(k,j,i) = (Fip*A1(k,j,i+1)-Fim*A1(k,j,i));
          divF(k,j,i) *= flux_pre;
          divF(k,j,i) /= dV(k,j,i);

        }
      }
    }
  } else if (kappatypeGlob=="usertable") {

    for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
          A1_out(k,j,i) = A1(k,j,i);

          real logtaum = std::log10(FMAX(variables["tau"](k,j,i-1)*idfx::units.GetDensity()*idfx::units.GetLength(),1.e-15));
          real Fim = pow(10.,kappatableGlob->GetHost(&logtaum))*A1(k,j,i)/std::pow(x1l(i),2.);
          real logtaup = std::log10(FMAX(variables["tau"](k,j,i)*idfx::units.GetDensity()*idfx::units.GetLength(),1.e-15));
          real Fip = pow(10.,kappatableGlob->GetHost(&logtaup))*A1(k,j,i+1)/std::pow(x1l(i+1),2.);
          divF(k,j,i)  = flux_pre*(Fip-Fim)/dV(k,j,i);

        }
      }
    }
  }

  Kokkos::deep_copy(variables["divF"], divF);
  Kokkos::deep_copy(variables["A1"], A1_out);

}

void InternalBoundary(Hydro *hydro, const real t) {
  auto *data = hydro->data;
  auto units = idfx::units;

  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray4D<real> Vs = hydro->Vs;
  IdefixArray1D<real> x1=data->x[IDIR];
  IdefixArray1D<real> x2=data->x[JDIR];

  real densityFloor = densityFloorGlob;
  idefix_for("InternalBoundary",
    0, data->np_tot[KDIR],
    0, data->np_tot[JDIR],
    0, data->np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                if(Vc(RHO,k,j,i)*units.GetDensity() <= densityFloor) {
                  real T= Vc(PRS,k,j,i)/Vc(RHO,k,j,i);
                  Vc(RHO,k,j,i)=densityFloor/units.GetDensity();
                  Vc(PRS,k,j,i)=T*Vc(RHO,k,j,i);
                  Vc(VX1,k,j,i)=ZERO_F;
                  Vc(VX2,k,j,i)=ZERO_F;
                }
              });
}

// Default constructor
// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output)
{
  auto Vc = data.hydro->Vc;

  // Mirror data on Host
  alphaMRIGlob = input.Get<real>("Setup","alphaMRI",0);
  alphaDZGlob = input.Get<real>("Setup","alphaDZ",0);
  epsilonGlob = input.Get<real>("Setup","epsilon",0);
  R0Glob = input.Get<real>("Setup","R0",0);
  rho0Glob = input.Get<real>("Setup","rho0",0);
  rhoindexGlob = input.Get<real>("Setup","rhoindex",0);
  rhominGlob = input.Get<real>("Setup","rhomin",0);
  T0Glob = input.Get<real>("Setup","T0",0);
  TsubGlob = input.Get<real>("Setup","Tsub",0);
  rhosGlob = input.Get<real>("Setup","rhos",0);
  rhosindexGlob = input.Get<real>("Setup","rhos_index",0);
  TwidthGlob = input.Get<real>("Setup","Twidth",0);
  TMriGlob = input.Get<real>("Setup","TMRI",0);
  f0Glob = input.Get<real>("Setup","f0",0);
  kappastarGlob = input.Get<real>("Setup","kappa_star",0);
  kappagasGlob = input.Get<real>("Setup","kappa_gas",0);
  densityFloorGlob = input.Get<real>("Setup","density_floor",0);
  alphawidthGlob = input.Get<real>("Setup","alphawidth",0);

  GGlob = input.Get<real>("Gravity","gravCst",0);
  MGlob = input.Get<real>("Gravity","Mcentral",0);

  muGlob = input.Get<real>("Hydro","mu",0);
  gammaGlob=data.hydro->eos->GetGamma();
  rsGlob=input.Get<real>("Rad","irr",1);
  TsGlob=input.Get<real>("Rad","irr",2);
  kappatypeGlob = input.Get<std::string>("Rad","irr",0);
  if (kappatypeGlob == "constant") {
    kappairrGlob = input.Get<real>("Rad","irr",3);
  } else if (kappatypeGlob == "usertable") {
    std::string irr_file = input.Get<std::string>("Rad","irr",4);
    kappatableGlob = new LookupTable<1>(irr_file,',');
  }

  columnGlob = new Column(IDIR,1,&data);
  columnGlob2 = new Column(IDIR,1,&data);

  auto temp_array = columnGlob->GetColumn();
  auto temp_array2 = columnGlob2->GetColumn();
  data.dump->RegisterVariable(temp_array,"tau");
  data.dump->RegisterVariable(temp_array2,"tau2");

  data.radiation[0]->EnrollKappa(&MyKappa);
  data.hydro->EnrollInternalBoundary(&InternalBoundary);
  //data.hydro->EnrollFluxBoundary(&FluxBoundary);
  data.hydro->viscosity->EnrollViscousDiffusivity(&MyViscosity);
  output.EnrollUserDefVariables(&ComputeUserVars);
  data.hydro->EnrollUserDefBoundary(&UserdefBoundaryNoStress);
  data.radiation[0]->EnrollUserDefBoundary(&UserdefBoundaryRad);

}

Setup::~Setup() {
  delete columnGlob;
}
