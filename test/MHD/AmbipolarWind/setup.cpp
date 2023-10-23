#include "idefix.hpp"
#include "setup.hpp"

real epsilonGlob;
real epsilonTopGlob;
real betaGlob;
real HidealGlob;
real AmMidGlob;
real gammaGlob;
real densityFloorGlob;
real trSmoothingGlob;


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

KOKKOS_INLINE_FUNCTION real computeVaMax(real t_change, real Va_ini_max, real Va_fin_max, real t){
  real Va_lim;

  if (t<t_change){
    Va_lim =t * (Va_fin_max - Va_ini_max)/t_change + Va_ini_max;
  }
  else{
    Va_lim = Va_fin_max;
  }

  return Va_lim;
}


void Ambipolar(DataBlock& data, real t, IdefixArray3D<real> &xAin) {
  IdefixArray3D<real> xA = xAin;
  IdefixArray1D<real> x1=data.x[IDIR];
  IdefixArray1D<real> x2=data.x[JDIR];
  IdefixArray4D<real> Vc=data.hydro->Vc;

  real Hideal = HidealGlob;
  real epsilon = epsilonGlob;
  real AmMid = AmMidGlob;
  real etamax = 10*epsilon*epsilon; // Corresponds to Rm=0.1
  real Rin = 1.0;
  real waveKillWidth = 0.1;
  real trSmoothing = trSmoothingGlob;

  idefix_for("Ambipolar",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real z=x1(i)*cos(x2(j));
                real R=FMAX(FABS(x1(i)*sin(x2(j))),ONE_F);
                real Omega=pow(R,-1.5);

                real zh = z/(R*epsilon);  // z in units of disc scale height h=R*epsilon
                real Am;

                  Am = AmMid / ( 0.5*(1-tanh((fabs(zh)-Hideal)/trSmoothing)));

                real B2 = Vc(BX1,k,j,i)*Vc(BX1,k,j,i)+Vc(BX2,k,j,i)*Vc(BX2,k,j,i)+Vc(BX3,k,j,i)*Vc(BX3,k,j,i);
                real eta = B2/(Omega*Am*Vc(RHO,k,j,i));
                if(eta>etamax) xA(k,j,i) = etamax/B2;
                else xA(k,j,i) = 1.0/(Omega*Am*Vc(RHO,k,j,i));

                // Kill it at the radial boundaryloop
                if(x1(i)/Rin < Rin*(1+waveKillWidth)) {
                  real w = (x1(i)-Rin)/(Rin*waveKillWidth);

                  xA(k,j,i) = xA(k,j,i)*w;
                }
              });

}

void Resistivity(DataBlock& data, real t, IdefixArray3D<real> &etain) {
  IdefixArray3D<real> eta = etain;
  IdefixArray1D<real> x1=data.x[IDIR];
  IdefixArray1D<real> x2=data.x[JDIR];
  IdefixArray4D<real> Vc=data.hydro->Vc;

  real epsilon = epsilonGlob;

  idefix_for("Resistivity",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                    eta(k,j,i) = x1(i) < 1.2 ? (1.2-x1(i))/0.2*epsilon*epsilon : ZERO_F;
              });

}

void MySourceTerm(Hydro *hydro, const real t, const real dtin) {
  auto *data = hydro->data;
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray4D<real> Uc = hydro->Uc;
  IdefixArray1D<real> x1=data->x[IDIR];
  IdefixArray1D<real> x2=data->x[JDIR];
  real epsilonTop = epsilonTopGlob;
  real epsilon = epsilonGlob;
  real tauGlob=0.1;
  real gamma_m1=gammaGlob-1.0;
  real dt=dtin;
  real Hideal=HidealGlob;
  real Rin=1.0;
  real trSmoothing = trSmoothingGlob;

  idefix_for("MySourceTerm",
    0, data->np_tot[KDIR],
    0, data->np_tot[JDIR],
    0, data->np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real r=x1(i);
                real th=x2(j);
                real z=r*cos(th);
                real R=r*sin(th);
                real R0=FMAX(R,Rin);
                real tau;

                real Zh = FABS(z/R)/epsilon;
                real Tdisk = epsilon*epsilon/R0;
                real Tcorona = epsilonTop*epsilonTop/R0;
                //if(x1(i) < 1.5) cscorona = csdisk;
                real Teff=0.5*(Tdisk+Tcorona)+0.5*(Tcorona-Tdisk)*tanh((Zh-Hideal)/trSmoothing);

                tau= tauGlob*(FMIN(pow(R,1.5),1.0));

                // Cooling /heatig function
                real Ptarget = Teff*Vc(RHO,k,j,i);

                Uc(ENG,k,j,i) += -dt*(Vc(PRS,k,j,i)-Ptarget)/(tau*gamma_m1);


                // inner shell relaxation
                /*
                if(r<1.2) {
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




void InternalBoundary(Hydro *hydro, const real t) {
  auto *data = hydro->data;
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray4D<real> Vs = hydro->Vs;
  IdefixArray1D<real> x1=data->x[IDIR];
  IdefixArray1D<real> x2=data->x[JDIR];

  real vAmax=computeVaMax(4.0,50.0,8.0,t);
  real densityFloor0 = densityFloorGlob;
  real Rin = 1.0;
  real epsilon=epsilonGlob;

  idefix_for("InternalBoundary",
    0, data->np_tot[KDIR],
    0, data->np_tot[JDIR],
    0, data->np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real R=x1(i)*sin(x2(j));
                real z=x1(i)*cos(x2(j));
                real zh = FABS(z/R)/epsilon;

                real b2=EXPAND(Vc(BX1,k,j,i)*Vc(BX1,k,j,i) , +Vc(BX2,k,j,i)*Vc(BX2,k,j,i), +Vc(BX3,k,j,i)*Vc(BX3,k,j,i) ) ;
                real va2=b2/Vc(RHO,k,j,i);
                real myMax=vAmax;
                //if(x1(i)<1.1) myMax=myMax/50.0;
                if(va2>myMax*myMax) {
                  real T = Vc(PRS,k,j,i)/Vc(RHO,k,j,i);
                  Vc(RHO,k,j,i) = b2/(myMax*myMax);
                  Vc(PRS,k,j,i) = T*Vc(RHO,k,j,i);
                }
                real densityFloor = computeDensityFloor(R,z,densityFloor0,Rin,epsilon);
                if(Vc(RHO,k,j,i) < densityFloor) {
                  real T= Vc(PRS,k,j,i)/Vc(RHO,k,j,i);
                  Vc(RHO,k,j,i)=densityFloor;
                  //Vc(PRS,k,j,i)=T*Vc(RHO,k,j,i);
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
void UserdefBoundary(Hydro *hydro, int dir, BoundarySide side, real t) {
    auto *data = hydro->data;
    if( (dir==IDIR) && (side == left)) {
        IdefixArray4D<real> Vc = hydro->Vc;
        IdefixArray4D<real> Vs = hydro->Vs;
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

                Vc(RHO,k,j,i) = 1.0/(Rin*sqrt(Rin))  * exp(1.0/ (csdisk*csdisk) * (1.0/sqrt(Rin*Rin+z*z)-1.0/Rin));
                real densityFloor = computeDensityFloor(R,z,densityFloor0,Rin,epsilon);
                if(Vc(RHO,k,j,i) < densityFloor) Vc(RHO,k,j,i) = densityFloor;

                Vc(PRS,k,j,i) = Vc(RHO,k,j,i)*csdisk*csdisk;

                if(Vc(VX1,k,j,ighost)>=ZERO_F) Vc(VX1,k,j,i) = -Vc(VX1,k,j,2*ighost-i);
                       else Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost);
                Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost);
                //real Rmin = FMAX(0.3,R);

                //Vc(VX3,k,j,i) = 1.0/sqrt(Rmin) * sqrt( Rmin / sqrt(Rmin*Rmin + z*z));
                Vc(VX3,k,j,i) = Omega*R;
                Vc(BX3,k,j,i) = - Vc(BX3,k,j,2*ighost-i);
                //Vc(BX3,k,j,i) = Vc(BX3,k,j,ighost);

            });
      hydro->boundary->BoundaryForX2s("UserDefX1",dir,side,
        KOKKOS_LAMBDA (int k, int j, int i) {
            Vs(BX2s,k,j,i) = Vs(BX2s,k,j,ighost);
          });

    }

    if( (dir==IDIR) && (side == right)) {
        IdefixArray4D<real> Vc = hydro->Vc;
        IdefixArray4D<real> Vs = hydro->Vs;
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
                Vc(BX3,k,j,i) = - Vc(BX3,k,j,2*ighost-i);
                //Vc(BX3,k,j,i) = Vc(BX3,k,j,ighost);

            });
      hydro->boundary->BoundaryForX2s("UserDefX1",dir,side,
        KOKKOS_LAMBDA (int k, int j, int i) {
            Vs(BX2s,k,j,i) = Vs(BX2s,k,j,ighost);
          });

    }


}

void EmfBoundary(DataBlock& data, const real t) {
    IdefixArray3D<real> Ex1 = data.hydro->emf->ex;
    IdefixArray3D<real> Ex2 = data.hydro->emf->ey;
    IdefixArray3D<real> Ex3 = data.hydro->emf->ez;
    if(data.lbound[IDIR] == userdef) {

        int ighost = data.beg[IDIR];

        idefix_for("EMFBoundary",0,data.np_tot[KDIR],0,data.np_tot[JDIR],
                    KOKKOS_LAMBDA (int k, int j) {
            Ex3(k,j,ighost) = ZERO_F;
        });
    }
    // additional zero EMF on the boundary
    if(data.lbound[JDIR] == axis) {
        int jghost = data.beg[JDIR];
        //printf("I'mbeing called\n");
        idefix_for("EMFBoundary",0,data.np_tot[KDIR],0,data.np_tot[IDIR],
                    KOKKOS_LAMBDA (int k, int i) {
            Ex3(k,jghost,i) = ZERO_F;
        });
    }
    if(data.rbound[JDIR] == axis) {
        int jghost = data.end[JDIR];
        //printf("I'mbeing called\n");
        idefix_for("EMFBoundary",0,data.np_tot[KDIR],0,data.np_tot[IDIR],
                    KOKKOS_LAMBDA (int k, int i) {
            Ex3(k,jghost,i) = ZERO_F;
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

void ComputeUserVars(DataBlock & data, UserDefVariablesContainer &variables) {

  // Use Invdt as scratch array
  IdefixArray3D<real> scrh("Scratch", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);

  // Ask for a computation of xA ambipolar in this scratch array
  Ambipolar(data, data.t, scrh);

  // Mirror data on Host
  DataBlockHost d(data);

  // Sync it
  d.SyncFromDevice();

  // Make references to the user-defined arrays (variables is a container of IdefixHostArray3D)
  // Note that the labels should match the variable names in the input file
  IdefixHostArray3D<real> Am  = variables["Am"];
  IdefixHostArray3D<real> InvDt  = variables["InvDt"];

  IdefixHostArray1D<real> x1=d.x[IDIR];
  IdefixHostArray1D<real> x2=d.x[JDIR];
  IdefixHostArray4D<real> Vc=d.Vc;
  IdefixArray3D<real>::HostMirror scrhHost = Kokkos::create_mirror_view(scrh);
  Kokkos::deep_copy(scrhHost,scrh);

  for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
    for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
      for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
        real z=x1(i)*cos(x2(j));
        real R=FMAX(FABS(x1(i)*sin(x2(j))),ONE_F);
        real Omega=pow(R,-1.5);
        Am(k,j,i) = 1.0/(Omega*scrhHost(k,j,i)*Vc(RHO,k,j,i));
        InvDt(k,j,i) = d.InvDt(k,j,i);
      }
    }
  }
}

// Default constructor


// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
    // Set the function for userdefboundary
    data.hydro->EnrollUserDefBoundary(&UserdefBoundary);
    data.hydro->EnrollAmbipolarDiffusivity(&Ambipolar);
    //data.hydro->EnrollOhmicDiffusivity(&Resistivity);
    data.hydro->EnrollUserSourceTerm(&MySourceTerm);
    data.hydro->EnrollInternalBoundary(&InternalBoundary);
    data.hydro->EnrollEmfBoundary(&EmfBoundary);
    //data.hydro->EnrollFluxBoundary(&FluxBoundary);
    output.EnrollUserDefVariables(&ComputeUserVars);
    gammaGlob=data.hydro->eos->GetGamma();
    epsilonGlob = input.Get<real>("Setup","epsilon",0);
    epsilonTopGlob = input.Get<real>("Setup","epsilonTop",0);
    betaGlob = input.Get<real>("Setup","beta",0);
    HidealGlob = input.Get<real>("Setup","Hideal",0);
    AmMidGlob = input.Get<real>("Setup","Am",0);
    densityFloorGlob = input.Get<real>("Setup","densityFloor",0);
    trSmoothingGlob = input.Get<real>("Setup","transitionSmoothing",0);
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


    real Rin=1.0;
    real m=-5.0/4.0;
    real B0 = epsilonGlob*sqrt(2.0/betaGlob);

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                real r=d.x[IDIR](i);
                real th=d.x[JDIR](j);
                real z=r*cos(th);
                real R=r*sin(th);
                if(R>Rin) {
                    real Zh = FABS(z/R)/epsilonGlob;
                    real csdisk = epsilonGlob/sqrt(R);
                    real cs2 = csdisk*csdisk;
                    d.Vc(RHO,k,j,i) = 1.0/(R*sqrt(R))  * exp(1.0/ (csdisk*csdisk) * (1.0/sqrt(R*R+z*z)-1.0/R)) ;
                    d.Vc(VX3,k,j,i) = 1.0/sqrt(R) * sqrt( FMAX(R / sqrt(R*R + z*z) -2.5*csdisk*csdisk,0.0) );
                    d.Vc(PRS,k,j,i) = cs2*d.Vc(RHO,k,j,i);
                    if(std::isnan(d.Vc(VX3,k,j,i))) {
                      idfx::cout << "Nan in R>Rin at (i,j,k)=("<< i <<"," << j << "," << k << "), (r,th,R,z)=(" << r << "," << th << "," << R << "," << z << ")" << std::endl;
                      IDEFIX_ERROR("Nan!s");
                    }
                }
                else {
                  real Zh = FABS(z/Rin)/epsilonGlob;
                  real csdisk = epsilonGlob/sqrt(Rin);
                  real cs2 = csdisk*csdisk;
                  d.Vc(RHO,k,j,i) = 1.0/(Rin*sqrt(Rin))  * exp(1.0/ (csdisk*csdisk) * (1.0/sqrt(Rin*Rin+z*z)-1.0/Rin));
                  d.Vc(VX3,k,j,i) = 1.0/sqrt(Rin) * sqrt( FMAX(Rin / sqrt(Rin*Rin + z*z) -2.5*csdisk*csdisk,0.0) );
                  d.Vc(PRS,k,j,i) = cs2*d.Vc(RHO,k,j,i);
                  if(std::isnan(d.Vc(VX3,k,j,i))) {
                    idfx::cout << "Nan in R<Rin at (i,j,k)=("<< i <<"," << j << "," << k << "), (r,th,R,z)=(" << r << "," << th << "," << R << "," << z << ")" << std::endl;
                    IDEFIX_ERROR("Nan!s");
                  }
                }

                d.Vc(VX1,k,j,i) = ZERO_F;
                d.Vc(VX2,k,j,i) = ZERO_F;

                real densityFloor = computeDensityFloor(R,z,densityFloorGlob,Rin,epsilonGlob);
                if(d.Vc(RHO,k,j,i) < densityFloor) {
                  d.Vc(RHO,k,j,i) = densityFloor;
                  //d.Vc(PRS,k,j,i) = T2*d.Vc(RHO,k,j,i);
                }

                // Vector potential on the corner
                real s=sin(d.xl[JDIR](j));
                R=d.xl[IDIR](i) * s;

                A(IDIR,k,j,i) = ZERO_F;
                A(JDIR,k,j,i) = ZERO_F;

                #ifdef EVOLVE_VECTOR_POTENTIAL
                  if(R>Rin) {
                    d.Ve(AX3e,k,j,i) = B0*(pow(Rin,m+2.0)/R * (-1.0/(m+2.0)) + pow(R,m+1.0)/(m+2.0) + Rin*Rin/(2.0*R));
                  }
                  else {
                    d.Ve(AX3e,k,j,i) = B0*R/2.0;
                  }
                #else
                  if(R>Rin) {
                    A(KDIR,k,j,i) = B0*(pow(Rin,m+2.0)/R * (-1.0/(m+2.0)) + pow(R,m+1.0)/(m+2.0));
                    A(KDIR,k,j,i) = B0*(pow(Rin,m+2.0)/R * (-1.0/(m+2.0)) + pow(R,m+1.0)/(m+2.0) + Rin*Rin/(2.0*R));
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
