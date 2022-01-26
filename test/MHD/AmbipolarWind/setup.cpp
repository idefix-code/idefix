#include "idefix.hpp"
#include "setup.hpp"
#include "boundaryloop.hpp"

real epsilonGlob;
real epsilonTopGlob;
real betaGlob;
real HidealGlob;
real AmMidGlob;
real gammaGlob;
real densityFloorGlob;

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


void Ambipolar(DataBlock& data, real t, IdefixArray3D<real> &xAin) {
  IdefixArray3D<real> xA = xAin;
  IdefixArray1D<real> x1=data.x[IDIR];
  IdefixArray1D<real> x2=data.x[JDIR];
  IdefixArray4D<real> Vc=data.hydro.Vc;

  real Hideal = HidealGlob;
  real epsilon = epsilonGlob;
  real AmMid = AmMidGlob;
  real etamax = 10*epsilon*epsilon; // Corresponds to Rm=0.1

  idefix_for("Ambipolar",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real z=x1(i)*cos(x2(j));
                real R=FMAX(FABS(x1(i)*sin(x2(j))),ONE_F);
                real Omega=pow(R,-1.5);

                real zh = z/(R*epsilon);  // z in units of disc scale height h=R*epsilon
                real Am = AmMid / ( 0.5*(1-tanh((fabs(zh)-Hideal)/0.5)));

                real B2 = Vc(BX1,k,j,i)*Vc(BX1,k,j,i)+Vc(BX2,k,j,i)*Vc(BX2,k,j,i)+Vc(BX3,k,j,i)*Vc(BX3,k,j,i);
                real eta = B2/(Omega*Am*Vc(RHO,k,j,i));
                if(eta>etamax) xA(k,j,i) = etamax/B2;
                else xA(k,j,i) = 1.0/(Omega*Am*Vc(RHO,k,j,i));
              });

}

void Resistivity(DataBlock& data, real t, IdefixArray3D<real> &etain) {
  IdefixArray3D<real> eta = etain;
  IdefixArray1D<real> x1=data.x[IDIR];
  IdefixArray1D<real> x2=data.x[JDIR];
  IdefixArray4D<real> Vc=data.hydro.Vc;

  real epsilon = epsilonGlob;

  idefix_for("Resistivity",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                    eta(k,j,i) = x1(i) < 1.2 ? (1.2-x1(i))/0.2*epsilon*epsilon : ZERO_F;
              });

}

void MySourceTerm(DataBlock &data, const real t, const real dtin) {
  IdefixArray4D<real> Vc = data.hydro.Vc;
  IdefixArray4D<real> Uc = data.hydro.Uc;
  IdefixArray1D<real> x1=data.x[IDIR];
  IdefixArray1D<real> x2=data.x[JDIR];
  real epsilonTop = epsilonTopGlob;
  real epsilon = epsilonGlob;
  real tauGlob=0.1;
  real gamma_m1=gammaGlob-1.0;
  real dt=dtin;
  real Hideal=HidealGlob;
  real tau2=0.5;
  idefix_for("MySourceTerm",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real r=x1(i);
                real th=x2(j);
                real z=r*cos(th);
                real R=r*sin(th);
                real cs2, tau;
                if(R>1.0) {
                    real Zh = FABS(z/R)/epsilon;
                    real csdisk = epsilon/sqrt(R);
                    real cscorona = epsilonTop/sqrt(R);
                    cs2=0.5*(csdisk*csdisk+cscorona*cscorona)+0.5*(cscorona*cscorona-csdisk*csdisk)*tanh(6*log(Zh/Hideal));
                    tau= tauGlob*sqrt(R);
                }
                else {
                    real Zh = FABS(z)/epsilon;
                    real csdisk = epsilon;
                    real cscorona = epsilonTop;
                    cs2=0.5*(csdisk*csdisk+cscorona*cscorona)+0.5*(cscorona*cscorona-csdisk*csdisk)*tanh(6*log(Zh/Hideal));
                    tau=tauGlob;
                }
                // Cooling /heatig function
                real Ptarget = cs2*Vc(RHO,k,j,i);

                Uc(ENG,k,j,i) += -dt*(Vc(PRS,k,j,i)-Ptarget)/(tau*gamma_m1);

        // Velocity relaxation
        if(r<1.5) {
            Uc(MX1,k,j,i) += -dt*(Vc(VX1,k,j,i)*Vc(RHO,k,j,i))/tau2;
            Uc(MX2,k,j,i) += -dt*(Vc(VX2,k,j,i)*Vc(RHO,k,j,i))/tau2;
        }

});


}



void InternalBoundary(DataBlock& data, const real t) {
  IdefixArray4D<real> Vc = data.hydro.Vc;
  IdefixArray4D<real> Vs = data.hydro.Vs;
  IdefixArray1D<real> x1=data.x[IDIR];
  IdefixArray1D<real> x2=data.x[JDIR];

  real vAmax=10.0;
  real densityFloor = densityFloorGlob;
  idefix_for("InternalBoundary",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real b2=EXPAND(Vc(BX1,k,j,i)*Vc(BX1,k,j,i) , +Vc(BX2,k,j,i)*Vc(BX2,k,j,i), +Vc(BX3,k,j,i)*Vc(BX3,k,j,i) ) ;
                real va2=b2/Vc(RHO,k,j,i);
                real myMax=vAmax;
                //if(x1(i)<1.1) myMax=myMax/50.0;
                if(va2>myMax*myMax) {
                  real T = Vc(PRS,k,j,i)/Vc(RHO,k,j,i);
                  Vc(RHO,k,j,i) = b2/(myMax*myMax);
                  Vc(PRS,k,j,i) = T*Vc(RHO,k,j,i);
                }
                if(Vc(RHO,k,j,i) < densityFloor) {
                  real T= Vc(PRS,k,j,i)/Vc(RHO,k,j,i);
                  Vc(RHO,k,j,i)=densityFloor;
                  Vc(PRS,k,j,i)=T*Vc(RHO,k,j,i);
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
void UserdefBoundary(DataBlock& data, int dir, BoundarySide side, real t) {

    if( (dir==IDIR) && (side == left)) {
        IdefixArray4D<real> Vc = data.hydro.Vc;
        IdefixArray4D<real> Vs = data.hydro.Vs;
        IdefixArray1D<real> x1 = data.x[IDIR];
        IdefixArray1D<real> x2 = data.x[JDIR];

        int ighost = data.nghost[IDIR];
        real Omega=1.0;
        data.hydro.boundary.BoundaryFor("UserDefX1",dir,side,
            KOKKOS_LAMBDA (int k, int j, int i) {
                real R=x1(i)*sin(x2(j));
                real z=x1(i)*cos(x2(j));

                Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost);
                Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost);
                if(Vc(VX1,k,j,ighost)>=ZERO_F) Vc(VX1,k,j,i) = - Vc(VX1,k,j,2*ighost-i);
                       else Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost);
                Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost);
                Vc(VX3,k,j,i) = 1.0*Omega; // 1.0=Rin
                Vs(BX2s,k,j,i) = Vs(BX2s,k,j,ighost);
                Vc(BX3,k,j,i) = ZERO_F;

            });
      data.hydro.boundary.BoundaryForX2s("UserDefX1",dir,side,
        KOKKOS_LAMBDA (int k, int j, int i) {
            Vs(BX2s,k,j,i) = Vs(BX2s,k,j,ighost);
          });

    }

}

void ComputeUserVars(DataBlock & data, UserDefVariablesContainer &variables) {

  // Use Invdt as scratch array
  IdefixArray3D<real> scrh = data.hydro.InvDt;

  // Ask for a computation of xA ambipolar in this scratch array
  Ambipolar(data, data.t, scrh);

  // Mirror data on Host
  DataBlockHost d(data);

  // Sync it
  d.SyncFromDevice();

  // Make references to the user-defined arrays (variables is a container of IdefixHostArray3D)
  // Note that the labels should match the variable names in the input file
  IdefixHostArray3D<real> Am  = variables["Am"];

  IdefixHostArray1D<real> x1=d.x[IDIR];
  IdefixHostArray1D<real> x2=d.x[JDIR];
  IdefixHostArray4D<real> Vc=d.Vc;
  IdefixHostArray3D<real> scrhHost=d.InvDt;

  for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
    for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
      for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
        real z=x1(i)*cos(x2(j));
        real R=FMAX(FABS(x1(i)*sin(x2(j))),ONE_F);
        real Omega=pow(R,-1.5);
        Am(k,j,i) = 1.0/(Omega*scrhHost(k,j,i)*Vc(RHO,k,j,i));
      }
    }
  }
}

// Default constructor


// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
    // Set the function for userdefboundary
    data.hydro.EnrollUserDefBoundary(&UserdefBoundary);
    data.hydro.EnrollAmbipolarDiffusivity(&Ambipolar);
    //data.hydro.EnrollOhmicDiffusivity(&Resistivity);
    data.hydro.EnrollUserSourceTerm(&MySourceTerm);
    data.hydro.EnrollInternalBoundary(&InternalBoundary);
    //data.hydro.EnrollEmfBoundary(&EmfBoundary);
    output.EnrollUserDefVariables(&ComputeUserVars);
    gammaGlob=data.hydro.GetGamma();
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
                    real cscorona = epsilonTopGlob/sqrt(R);
                    real cs2=0.5*(csdisk*csdisk+cscorona*cscorona)+0.5*(cscorona*cscorona-csdisk*csdisk)*tanh(6*log(Zh/HidealGlob));
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
                  real cscorona = epsilonTopGlob/sqrt(Rin);
                  real cs2=0.5*(csdisk*csdisk+cscorona*cscorona)+0.5*(cscorona*cscorona-csdisk*csdisk)*tanh(6*log(Zh/HidealGlob));
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

                if(d.Vc(RHO,k,j,i) < densityFloorGlob) {
                  real T2=d.Vc(PRS,k,j,i)/d.Vc(RHO,k,j,i);
                  d.Vc(RHO,k,j,i) = densityFloorGlob;
                  d.Vc(PRS,k,j,i) = T2*d.Vc(RHO,k,j,i);
                }

                // Vector potential on the corner
                real s=sin(d.xl[JDIR](j));
                R=d.xl[IDIR](i) * s;

                A(IDIR,k,j,i) = ZERO_F;
                A(JDIR,k,j,i) = ZERO_F;

                if(R>Rin) {
                  A(KDIR,k,j,i) = B0*(pow(Rin,m+2.0)/R * (-1.0/(m+2.0)) + pow(R,m+1.0)/(m+2.0));
                }
                else {
                  A(KDIR,k,j,i) = 0.0;
                }

            }
        }
    }

    // Make the field from the vector potential
    d.MakeVsFromAmag(A);


    // Send it all, if needed
    d.SyncToDevice();
}
