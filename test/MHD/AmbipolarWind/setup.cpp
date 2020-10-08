#include "idefix.hpp"
#include "setup.hpp"

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
  IdefixArray4D<real> Vc=data.Vc;

  real Hideal = HidealGlob;
  real epsilon = epsilonGlob;
  real AmMid = AmMidGlob;
  real etamax = 10*epsilon*epsilon; // Corresponds to Rm=0.1

  idefix_for("Ambipolar",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real z=x1(i)*cos(x2(j));
                real R=FMAX(FABS(x1(i)*sin(x2(j))),ONE_F);
                real Omega=pow(R,-1.5);

                real q = z/(Hideal*R*epsilon);
                real Am = AmMid * exp(q*q*q*q);

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
  IdefixArray4D<real> Vc=data.Vc;

  real epsilon = epsilonGlob;

  idefix_for("Resistivity",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                    eta(k,j,i) = x1(i) < 1.2 ? epsilon*epsilon : ZERO_F;
              });

}

void MySourceTerm(DataBlock &data, const real t, const real dtin) {
  IdefixArray4D<real> Vc = data.Vc;
  IdefixArray4D<real> Uc = data.Uc;
  IdefixArray1D<real> x1=data.x[IDIR];
  IdefixArray1D<real> x2=data.x[JDIR];
  real epsilonTop = epsilonTopGlob;
  real epsilon = epsilonGlob;
  real tauGlob=0.1;
  real gamma_m1=gammaGlob-1.0;
  real dt=dtin;
  real Hideal=HidealGlob;
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
		if(R<1.2) {
			Uc(MX1,k,j,i) += -dt*(Vc(VX1,k,j,i)*Vc(RHO,k,j,i));
			Uc(MX2,k,j,i) += -dt*(Vc(VX2,k,j,i)*Vc(RHO,k,j,i));
		}

});


}

void EmfBoundary(DataBlock& data, const real t) {
    IdefixArray3D<real> Ex1 = data.emf.ex;
    IdefixArray3D<real> Ex2 = data.emf.ey;
    IdefixArray3D<real> Ex3 = data.emf.ez;
    if(data.lbound[IDIR] == userdef) {

        int ighost = data.nghost[IDIR];

        idefix_for("EMFBoundary",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,ighost+1,
                    KOKKOS_LAMBDA (int k, int j, int i) {
            Ex3(k,j,i) = ZERO_F;
        });
    }
    if(data.lbound[JDIR] == userdef) {
        int jghost = data.nghost[JDIR];
        //printf("I'mbeing called\n");
        idefix_for("EMFBoundary",0,data.np_tot[KDIR],0,data.np_tot[IDIR],
                    KOKKOS_LAMBDA (int k, int i) {
            Ex3(k,jghost,i) = ZERO_F;
        });
    }
    if(data.rbound[JDIR] == userdef) {
        int jghost = data.end[JDIR];
        //printf("I'mbeing called\n");
        idefix_for("EMFBoundary",0,data.np_tot[KDIR],0,data.np_tot[IDIR],
                    KOKKOS_LAMBDA (int k, int i) {
            Ex3(k,jghost,i) = ZERO_F;
        });
    }
}

void InternalBoundary(DataBlock& data, const real t) {
  IdefixArray4D<real> Vc = data.Vc;
  IdefixArray4D<real> Vs = data.Vs;
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
        IdefixArray4D<real> Vc = data.Vc;
        IdefixArray4D<real> Vs = data.Vs;
        IdefixArray1D<real> x1 = data.x[IDIR];
        IdefixArray1D<real> x2 = data.x[JDIR];

        int ighost = data.nghost[IDIR];
        real Omega=1.0;
        idefix_for("UserDefBoundary",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,ighost,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        real R=x1(i)*sin(x2(j));
                        real z=x1(i)*cos(x2(j));

                        Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost);
                        Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost);
                        if(Vc(VX1,k,j,ighost)>=ZERO_F) Vc(VX1,k,j,i) = - Vc(VX1,k,j,2*ighost-i);
			                   else Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost);
                        Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost);
                        Vc(VX3,k,j,i) = R*Omega;
                        Vs(BX2s,k,j,i) = Vs(BX2s,k,j,ighost);
                        Vc(BX3,k,j,i) = ZERO_F;

                    });

    }

    if( dir==JDIR) {
        IdefixArray4D<real> Vc = data.Vc;
        IdefixArray4D<real> Vs = data.Vs;
        int jghost;
        int jbeg,jend;
        if(side == left) {
            jghost = data.beg[JDIR];
            jbeg = 0;
            jend = data.beg[JDIR];
            //return;
        }
        else if(side==right) {
            jghost = data.end[JDIR]-1;
            jbeg=data.end[JDIR];
            jend=data.np_tot[JDIR];
        }


        idefix_for("UserDefBoundary",0,data.np_tot[KDIR],jbeg,jend,0,data.np_tot[IDIR],
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        Vc(RHO,k,j,i) = Vc(RHO,k,jghost,i);
                        Vc(PRS,k,j,i) = Vc(PRS,k,jghost,i);
                        Vc(VX1,k,j,i) = Vc(VX1,k,jghost,i);
                        Vc(VX2,k,j,i) = - Vc(VX2,k,2*jghost-j,i);
                        Vc(VX3,k,j,i) = ZERO_F;
                        Vs(BX1s,k,j,i) = Vs(BX1s,k,jghost,i);
                        Vc(BX3,k,j,i) = ZERO_F;

                    });


    }

}

void Potential(DataBlock& data, const real t, IdefixArray1D<real>& x1, IdefixArray1D<real>& x2, IdefixArray1D<real>& x3, IdefixArray3D<real>& phi) {

    idefix_for("Potential",0,data.np_tot[KDIR], 0, data.np_tot[JDIR], 0, data.np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
        phi(k,j,i) = -1.0/x1(i);
    });

}


// Default constructor
Setup::Setup() {}

// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Hydro &hydro) {
    // Set the function for userdefboundary
    hydro.EnrollUserDefBoundary(&UserdefBoundary);
    hydro.EnrollGravPotential(&Potential);
    hydro.EnrollAmbipolarDiffusivity(&Ambipolar);
    //hydro.EnrollOhmicDiffusivity(&Resistivity);
    hydro.EnrollUserSourceTerm(&MySourceTerm);
    hydro.EnrollInternalBoundary(&InternalBoundary);
    hydro.EnrollEmfBoundary(&EmfBoundary);
    hydro.SetGamma(1.05);
    gammaGlob=hydro.GetGamma();
    epsilonGlob = input.GetReal("Setup","epsilon",0);
    epsilonTopGlob = input.GetReal("Setup","epsilonTop",0);
    betaGlob = input.GetReal("Setup","beta",0);
    HidealGlob = input.GetReal("Setup","Hideal",0);
    AmMidGlob = input.GetReal("Setup","Am",0);
    densityFloorGlob = input.GetReal("Setup","densityFloor",0);
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
                    d.Vc(VX3,k,j,i) = 1.0/sqrt(R) * sqrt( R / sqrt(R*R + z*z) -2.5*csdisk*csdisk );
                    d.Vc(PRS,k,j,i) = cs2*d.Vc(RHO,k,j,i);

                }
                else {
                  real Zh = FABS(z/Rin)/epsilonGlob;
                  real csdisk = epsilonGlob/sqrt(Rin);
                  real cscorona = epsilonTopGlob/sqrt(Rin);
                  real cs2=0.5*(csdisk*csdisk+cscorona*cscorona)+0.5*(cscorona*cscorona-csdisk*csdisk)*tanh(6*log(Zh/HidealGlob));
                  d.Vc(RHO,k,j,i) = 1.0/(Rin*sqrt(Rin))  * exp(1.0/ (csdisk*csdisk) * (1.0/sqrt(Rin*Rin+z*z)-1.0/Rin));
                  d.Vc(VX3,k,j,i) = 1.0/sqrt(Rin) * sqrt( Rin / sqrt(Rin*Rin + z*z) -2.5*csdisk*csdisk );
                  d.Vc(PRS,k,j,i) = cs2*d.Vc(RHO,k,j,i);
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

// Analyse data to produce an output
void Setup::MakeAnalysis(DataBlock & data, real t) {

}




// Do a specifically designed user step in the middle of the integration
void ComputeUserStep(DataBlock &data, real t, real dt) {

}
