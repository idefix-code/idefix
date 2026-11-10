#include "idefix.hpp"
#include "setup.hpp"

real RICGlob;
real RinGlob;
real rho0Glob;
real T0Glob;
real TXGlob;
real LstarGlob;
real kappaXGlob;
real muGlob;
real gammaGlob;
real MGlob;
real GGlob;
real densityFloorGlob;
real rhoindexGlob;
real x1begGlob;
real epsilonGlob;
real kappapGlob;
real kapparGlob;
real kramersTindexGlob;
real kramersrhoindexGlob;
real kramerspconstantGlob;
real kramersrconstantGlob;
real tstartGlob;
real ttransGlob;

std::string kappatypeGlob;

LookupTable<2> *kappaptabGlob;
LookupTable<2> *kappartabGlob;

bool haveRadiationGlob = false;

Column *column;

void MySourceTerm(Fluid<DefaultPhysics> *hydro, const real tin, const real dtin) {
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray4D<real> Uc = hydro->Uc;
  auto data = hydro->data;
  auto units = idfx::units;
  IdefixArray3D<real> InvDt = hydro->InvDt;
  IdefixArray1D<real> r = data->x[IDIR];
  real mu = muGlob;
  real TX = TXGlob;
  real kappaX = kappaXGlob;
  real Rin = RinGlob;
  real Lstar = LstarGlob;

  real tstart = tstartGlob;
  real t_trans = ttransGlob;
  real dt = dtin;
  real t = tin ;

  IdefixArray3D<real> tau;
  IdefixArray3D<real> kapparho("rho",data->np_tot[KDIR],data->np_tot[JDIR],data->np_tot[IDIR]);
  idefix_for("Compute kapparho in MySourceTerm",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA(int k, int j, int i) {
      kapparho(k,j,i) = Vc(RHO,k,j,i)*kappaX*units.GetDensity()*units.GetLength();
    });
  column->ComputeColumn(kapparho);
  tau = column->GetColumn();

  idefix_for("MySourceTerm",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {

                real FX = Lstar * std::exp(-tau(k,j,i));
                FX /= r(i)*r(i)*units.GetLength()*units.GetLength();
                real nH = Vc(RHO,k,j,i)*units.GetDensity()/(units.u*mu);

                real xi = FX/nH;

                real T = Vc(PRS,k,j,i)/Vc(RHO,k,j,i)*units.GetKelvin()*mu;

                real ne1 = 2.545825721359034e-52*std::pow(T,12.27740153);
                real ne2 = 0.00015577691900385798*std::pow(T,0.86092628);
                real ne3 = 1.21;
                real ne = FMIN(FMIN(ne1,ne2),ne3);
                ne *= nH;

                //Blondin formula from Proga, Stone & Kallman 2000
                //Compton heating/cooling
                real Gcompton = ne*8.9e-36*xi*(TX-4.*T); //Compton heating/cooling
                //if (std::isnan(Gcompton)) {
                //  std::printf("xi=%e, T=%e\n",xi,T);
                //  Kokkos::abort("Nan in Gcompton");
                //}
                //net rate of X-ray photoionization heating-recombination cooling
                real GX = nH*1.5e-21*std::pow(xi,0.25)*std::pow(T,-0.5)*(1.-T/TX);
                //if (std::isnan(GX)) {
                //  std::printf("xi=%e, T=%e, PRS=%e, RHO=%e\n",xi,T,Vc(PRS,k,j,i),Vc(RHO,k,j,i));
                //  Kokkos::abort("Nan in GX");
                //}
                //Bremsstrahlung cooling
                real Lbt = ne*3.3e-27*std::sqrt(T);
                //Line cooling
                real delta = 0.; //parameter controlling line cooling
        delta *= FMIN(ONE_F,(t*units.GetTime()-tstart)/t_trans); // don't start line cooling at beginning
                Lbt += ne*(1.7e-18*std::exp(FMAX(-1.3e5/T,-30.))/xi*std::pow(T,-0.5)+1.e-24)*delta;

                real total_cooling = nH*(Gcompton+GX-Lbt)*units.GetTime()/units.GetEnergy();
                total_cooling *= FMIN(ONE_F,(t*units.GetTime()-tstart)/t_trans);

                real new_dt = dt;
                real frac_ENG = 0.8;
                if (total_cooling*dt > frac_ENG*Uc(ENG,k,j,i)) {
                    InvDt(k,j,i) = total_cooling/(frac_ENG*Uc(ENG,k,j,i));
                    new_dt = 1./InvDt(k,j,i);
                }
                Uc(ENG,k,j,i) -= new_dt*total_cooling;
              }
            );
}


void InternalBoundary(Hydro *hydro, const real t) {
  auto *data = hydro->data;
  auto units = idfx::units;

  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray4D<real> Uc = hydro->Uc;
  IdefixArray4D<real> Vs = hydro->Vs;
  IdefixArray1D<real> x1=data->x[IDIR];
  IdefixArray1D<real> x2=data->x[JDIR];

  real densityFloor = densityFloorGlob;
  real rhoindex = rhoindexGlob;
  real x1beg = x1begGlob;
  real mu = muGlob;
  real gamma = gammaGlob;

  idefix_for("InternalBoundary",
    0, data->np_tot[KDIR],
    0, data->np_tot[JDIR],
    0, data->np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real densfloor = densityFloor*std::pow(x1beg/x1(i),rhoindex);
                if(Vc(RHO,k,j,i)*units.GetDensity() <= densfloor) {
                  real T= Vc(PRS,k,j,i)/Vc(RHO,k,j,i);
                  Vc(RHO,k,j,i)=densfloor/units.GetDensity();
                  Vc(PRS,k,j,i)=T*Vc(RHO,k,j,i);
                  Vc(VX1,k,j,i)=ZERO_F;
                  Vc(VX2,k,j,i)=ZERO_F;
                }
              });
}


void UserdefBoundary(Fluid<DefaultPhysics> *hydro, int dir, BoundarySide side, real t) {
  IdefixArray4D<real> Vc = hydro->Vc;
  auto *data = hydro->data;
  auto units = idfx::units;

  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];

  real rho0 = rho0Glob;
  real Rin = RinGlob;
  real RIC = RICGlob;
  real rhoindex = rhoindexGlob;
  real CG = MGlob*GGlob;
  real epsilon = epsilonGlob;
  real T0 = T0Glob;
  real mu = muGlob;

  if(dir==JDIR) {
    int jend,jbeg;
    if (side ==right) {
      jbeg = data->end[JDIR];
      jend =data->np_tot[JDIR];
      idefix_for("UserDefBoundary",
        0, data->np_tot[KDIR],
        jbeg, jend,
        0, data->np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {

          real R = FMAX(x1(i)*std::sin(x2(j)),Rin);
          real Omega = std::sqrt(CG)*std::pow(R,-1.5);
          //real vk = Omega*R;
          //real cs = epsilon*vk;
          real T = T0;
          real cs = std::sqrt(T/(units.GetKelvin()*mu));;

          Vc(RHO,k,j,i) = rho0*std::pow(RIC/R,rhoindex)/units.GetDensity();
          Vc(VX1,k,j,i) = ZERO_F;
          Vc(VX2,k,j,i) = ZERO_F;
          Vc(VX3,k,j,i) =  Omega*R;
          Vc(PRS,k,j,i) = Vc(RHO,k,j,i)*cs*cs;

        });
    }
  }
}


void FluxBoundary(Fluid<DefaultPhysics> *hydro, int dir, BoundarySide side, const real t) {
    idfx::pushRegion("FluxInternal");

    auto *data = hydro->data;

    IdefixArray4D<real> Flux = data->hydro->FluxRiemann[dir];
    IdefixArray4D<real> Vc = hydro->Vc;

    if( dir== IDIR ) {
      int ibeg,iend;
      if ( side == left ) {
        ibeg =0;
        iend = data->beg[IDIR];
        idefix_for("FluxIDIR",
                    0, data->np_tot[KDIR],
                    0, data->np_tot[JDIR],
                    ibeg, iend,
           KOKKOS_LAMBDA (int k, int j, int i) {
           if(Vc(VX1,k,j,i) >= ZERO_F) {
             Flux(RHO, k, j, i) = ZERO_F;
             Flux(MX1, k, j, i) = ZERO_F;
             Flux(MX2, k, j, i) = ZERO_F;
             Flux(MX3, k, j, i) = ZERO_F;
             Flux(ENG, k, j, i) = ZERO_F;
           }
        });
      }
    } else if ( dir == JDIR ) {
      int jbeg,jend;
      if ( side == left ) {
        jbeg = 0;
        jend = data->end[JDIR];
        idefix_for("FluxJDIR",
                    0, data->np_tot[KDIR],
                    jbeg, jend,
                    0, data->np_tot[IDIR],
           KOKKOS_LAMBDA (int k, int j, int i) {
           if(Vc(VX2,k,j,i) >= ZERO_F) {
             Flux(RHO, k, j, i) = ZERO_F;
             Flux(MX1, k, j, i) = ZERO_F;
             Flux(MX2, k, j, i) = ZERO_F;
             Flux(MX3, k, j, i) = ZERO_F;
             Flux(ENG, k, j, i) = ZERO_F;
           }
        });
      }
    }



    idfx::popRegion();

}

// Compute user variables which will be written in vtk files
void ComputeUserVars(DataBlock & data, UserDefVariablesContainer &variables) {
  // Mirror data on Host
  DataBlockHost d(data);

  // Sync it
  d.SyncFromDevice();

  auto units = idfx::units;

  IdefixArray4D<real> Vc=(&data)->hydro->Vc;

  IdefixHostArray3D<real> FX  = variables["FX"];
  IdefixHostArray3D<real> xi  = variables["xi"];
  IdefixHostArray3D<real> Hcompton  = variables["Hcompton"];
  IdefixHostArray3D<real> Ccompton  = variables["Ccompton"];
  IdefixHostArray3D<real> HX  = variables["HX"];
  IdefixHostArray3D<real> CX  = variables["CX"];
  IdefixHostArray3D<real> bremLbt  = variables["bremLbt"];
  IdefixHostArray3D<real> lineLbt  = variables["lineLbt"];
  IdefixHostArray3D<real> ne  = variables["ne"];

  real mu = muGlob;
  real TX = TXGlob;
  real kappaX = kappaXGlob;
  real Rin = RinGlob;
  real Lstar = LstarGlob;

  IdefixArray3D<real> tau;
  IdefixArray3D<real> kapparho("rho",data.np_tot[KDIR],data.np_tot[JDIR],data.np_tot[IDIR]);
  idefix_for("Compute kapparho in MySourceTerm",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
    KOKKOS_LAMBDA(int k, int j, int i) {
      kapparho(k,j,i) = Vc(RHO,k,j,i)*kappaX*units.GetDensity()*units.GetLength();
    });
  column->ComputeColumn(kapparho);
  tau = column->GetColumn();

  IdefixArray3D<real>::host_mirror_type tauHost;
  tauHost = Kokkos::create_mirror_view(tau);
  Kokkos::deep_copy(tauHost,tau);

  for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
    for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
      for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {

          FX(k,j,i) = Lstar * std::exp(-tauHost(k,j,i));
          FX(k,j,i) /= d.x[IDIR](i)*d.x[IDIR](i)*units.GetLength()*units.GetLength();
          real nH = d.Vc(RHO,k,j,i)*units.GetDensity()/(units.u*mu);
          xi(k,j,i) = FX(k,j,i)/nH;

          real T = d.Vc(PRS,k,j,i)/d.Vc(RHO,k,j,i)*units.GetKelvin()*mu;

          real ne1 = 2.545825721359034e-52*std::pow(T,12.27740153);
          real ne2 = 0.00015577691900385798*std::pow(T,0.86092628);
          real ne3 = 1.21;
          ne(k,j,i) = FMIN(FMIN(ne1,ne2),ne3);
          ne(k,j,i) *= nH;

          //Blondin formula from Proga, Stone & Kallman 2000
          //Compton heating/cooling
          Hcompton(k,j,i) = ne(k,j,i)*8.9e-36*xi(k,j,i)*TX; //Compton heating/cooling
          Ccompton(k,j,i) = ne(k,j,i)*8.9e-36*xi(k,j,i)*4.*T; //Compton heating/cooling
          //net rate of X-ray photoionization heating-recombination cooling
          HX(k,j,i) = nH*1.5e-21*std::pow(xi(k,j,i),0.25)*std::pow(T,-0.5);
          CX(k,j,i) = nH*1.5e-21*std::pow(xi(k,j,i),0.25)*std::pow(T,-0.5)*T/TX;
          //Bremsstrahlung cooling
          bremLbt(k,j,i) = ne(k,j,i)*3.3e-27*std::sqrt(T);
          //Line cooling
          real delta = 1.; //parameter controlling line cooling
          lineLbt(k,j,i) = ne(k,j,i)*(1.7e-18*std::exp(FMIN(-1.3e5/T,-30.))/xi(k,j,i)*std::pow(T,-0.5)+1.e-24)*delta;
      }
    }
  }

  Kokkos::deep_copy(variables["ne"],ne);
  Kokkos::deep_copy(variables["tau"],tau);
  Kokkos::deep_copy(variables["FX"], FX);
  Kokkos::deep_copy(variables["xi"], xi);
  Kokkos::deep_copy(variables["Hcompton"], Hcompton);
  Kokkos::deep_copy(variables["Ccompton"], Ccompton);
  Kokkos::deep_copy(variables["HX"], HX);
  Kokkos::deep_copy(variables["CX"], CX);
  Kokkos::deep_copy(variables["bremLbt"], bremLbt);
  Kokkos::deep_copy(variables["lineLbt"], lineLbt);

}

// Default constructor
// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output)
{
  auto Vc = data.hydro->Vc;

  // Mirror data on Host
  RICGlob = input.Get<real>("Setup","RIC",0);
  RinGlob = input.Get<real>("Grid","X1-grid",1);
  rho0Glob = input.Get<real>("Setup","rho0",0);
  rhoindexGlob = input.Get<real>("Setup","rhoindex",0);
  T0Glob = input.Get<real>("Setup","T0",0);
  densityFloorGlob = input.Get<real>("Setup","density_floor",0);
  epsilonGlob = input.Get<real>("Setup","epsilon",0);
  TXGlob = input.Get<real>("Setup","TX",0);
  LstarGlob = input.Get<real>("Setup","Lstar",0);
  kappaXGlob = input.Get<real>("Setup","kappaX",0);
  tstartGlob = input.Get<real>("Setup","tstart",0);
  ttransGlob = input.Get<real>("Setup","ttrans",0);

  GGlob = input.Get<real>("Gravity","gravCst",0);
  MGlob = input.Get<real>("Gravity","Mcentral",0);

  muGlob = input.Get<real>("Hydro","mu",0);
  gammaGlob=data.hydro->eos->GetGamma();

  x1begGlob = grid.xbeg[IDIR];

  column = new Column(IDIR,1,&data);

  data.hydro->EnrollInternalBoundary(&InternalBoundary);
  //data.hydro->EnrollFluxBoundary(&FluxBoundary);
  data.hydro->EnrollUserSourceTerm(&MySourceTerm);
  data.hydro->EnrollUserDefBoundary(&UserdefBoundary);
  output.EnrollUserDefVariables(&ComputeUserVars);

  if (input.CheckBlock("Rad")) {
    haveRadiationGlob = true;
    kappatypeGlob = input.Get<std::string>("Rad","kappa",0);
    if (kappatypeGlob.compare("usertable") == 0){
      std::string kappapfile = input.Get<std::string>("Rad","kappa",2);
      std::string kapparfile = input.Get<std::string>("Rad","kappa",3);
      kappaptabGlob = new LookupTable<2>(kappapfile,',');
      kappartabGlob = new LookupTable<2>(kapparfile,',');
    } else if (kappatypeGlob.compare("kramers") == 0){
      kramersTindexGlob = input.Get<real>("Rad","kappa",4);
      kramersrhoindexGlob = input.Get<real>("Rad","kappa",3);
      kramerspconstantGlob = input.Get<real>("Rad","kappa",1);
      kramersrconstantGlob = input.Get<real>("Rad","kappa",2);
    } else if (kappatypeGlob.compare("constant") == 0) {
      kappapGlob = input.Get<real>("Rad","kappa",1);
      kapparGlob = input.Get<real>("Rad","kappa",2);
    }
  }
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);

    auto units=idfx::units;

    real Rin = RinGlob;
    real RIC = RICGlob;
    real rho0 = rho0Glob;
    real densityFloor = densityFloorGlob;
    real rhoindex = rhoindexGlob;
    real epsilon = epsilonGlob;
    real T0 = T0Glob;
    real mu = muGlob;
    real CG = MGlob*GGlob;
    real x1beg = x1begGlob;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {

              real R = FMAX(d.x[IDIR](i)*std::sin(d.x[JDIR](j)),Rin);
              real z2 = std::pow(d.x[IDIR](i)*std::cos(d.x[JDIR](j)),2.);
              real Omega = std::sqrt(CG)*std::pow(R,-1.5);
              real vk = Omega*R;

              real T = T0;
              //real cs = epsilon*vk;
              real cs = std::sqrt(T/(units.GetKelvin()*mu));;
              real H = cs/Omega;

              real densfloor = densityFloor*std::pow(x1beg/d.x[IDIR](i),rhoindex);
              d.Vc(RHO,k,j,i) = (rho0*std::pow(RIC/R,rhoindex)*std::exp(-0.25*M_PI*z2/(H*H))+densfloor)/idfx::units.GetDensity();
              d.Vc(PRS,k,j,i) = d.Vc(RHO,k,j,i)*cs*cs;
              d.Vc(VX1,k,j,i) = 0.;
              d.Vc(VX2,k,j,i) = 0.;
              d.Vc(VX3,k,j,i) =  Omega*R;

              //real T = d.Vc(PRS,k,j,i)/d.Vc(RHO,k,j,i)*units.GetKelvin()*mu;
              if (haveRadiationGlob) {
                d.RadVc[0](ER,k,j,i) = idfx::units.ar*std::pow(T,4.)/idfx::units.GetEnergy();
                d.RadVc[0](FR1,k,j,i) = ZERO_F;
                d.RadVc[0](FR2,k,j,i) = ZERO_F;
                d.RadVc[0](FR3,k,j,i) = ZERO_F;
              }
            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}
