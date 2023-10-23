#include "idefix.hpp"
#include "setup.hpp"


real epsilonGlob;
real betaGlob;
real HidealGlob;
real AmMidGlob;
real gammaGlob;
real densityFloorGlob;


void MySourceTerm(Hydro *hydro, const real t, const real dtin) {
  auto *data = hydro->data;
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray4D<real> Uc = hydro->Uc;
  IdefixArray1D<real> x1=data->x[IDIR];
  IdefixArray1D<real> x2=data->x[JDIR];
  real epsilon = epsilonGlob;
  real tauGlob=0.1;
  real gamma_m1=gammaGlob-1.0;
  real dt=dtin;
  idefix_for("MySourceTerm",
    0, data->np_tot[KDIR],
    0, data->np_tot[JDIR],
    0, data->np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real r=x1(i);
                real th=x2(j);
                real z=r*cos(th);
                real R=r*sin(th);
                real cs2, tau;
                if(R>1.0) {
                    real csdisk = epsilon/sqrt(R);
                    cs2=(csdisk*csdisk);
                    tau= tauGlob*sqrt(R);
                }
                else {
                    real csdisk = epsilon;
                    cs2=(csdisk*csdisk);
                    tau=tauGlob;
                }
                // Cooling /heatig function
                real Ptarget = cs2*Vc(RHO,k,j,i);

                Uc(ENG,k,j,i) += -dt*(Vc(PRS,k,j,i)-Ptarget)/(tau*gamma_m1);

        // Velocity relaxation
            });


}

void FargoVelocity(DataBlock &data, IdefixArray2D<real> &Vphi) {
  IdefixArray1D<real> x1 = data.x[IDIR];
  IdefixArray1D<real> x2 = data.x[JDIR];

  idefix_for("FargoVphi",0,data.np_tot[JDIR], 0, data.np_tot[IDIR],
      KOKKOS_LAMBDA (int j, int i) {
      Vphi(j,i) = 1.0/sqrt(x1(i)*sin(x2(j)));
  });
}


void InternalBoundary(Hydro *hydro, const real t) {
  auto *data = hydro->data;
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray4D<real> Vs = hydro->Vs;
  IdefixArray1D<real> x1=data->x[IDIR];
  IdefixArray1D<real> x2=data->x[JDIR];

  real vAmax=10.0;
  real densityFloor = densityFloorGlob;
  idefix_for("InternalBoundary",
    0, data->np_tot[KDIR],
    0, data->np_tot[JDIR],
    0, data->np_tot[IDIR],
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
void UserdefBoundary(Hydro *hydro, int dir, BoundarySide side, real t) {
    auto *data = hydro->data;
    if( (dir==IDIR) && (side == left)) {
        IdefixArray4D<real> Vc = hydro->Vc;
        IdefixArray4D<real> Vs = hydro->Vs;
        IdefixArray1D<real> x1 = data->x[IDIR];
        IdefixArray1D<real> x2 = data->x[JDIR];

        int ighost = data->nghost[IDIR];
        real Omega=1.0;
        idefix_for("UserDefBoundary",
          0, data->np_tot[KDIR],
          0, data->np_tot[JDIR],
          0, ighost,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        real R=x1(i)*sin(x2(j));
                        real z=x1(i)*cos(x2(j));

                        Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost);
                        Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost);
                        if(Vc(VX1,k,j,ighost)>=ZERO_F) Vc(VX1,k,j,i) = - Vc(VX1,k,j,2*ighost-i);
                               else Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost);
                        Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost);
                        Vc(VX3,k,j,i) = R*Omega;
                    });

        idefix_for("UserDefBoundaryX1S",
          0, data->np_tot[KDIR],
          0, data->np_tot[JDIR]+1,
          0, ighost,
                      KOKKOS_LAMBDA (int k, int j, int i) {

                          Vs(BX2s,k,j,i) = ZERO_F;


                      });

        idefix_for("UserDefBoundaryX3S",
          0, data->np_tot[KDIR]+1,
          0, data->np_tot[JDIR],
          0, ighost,
                      KOKKOS_LAMBDA (int k, int j, int i) {
                          Vs(BX3s,k,j,i) = ZERO_F;

            });

    }

    if( dir==JDIR) {
        IdefixArray4D<real> Vc = hydro->Vc;
        IdefixArray4D<real> Vs = hydro->Vs;
        int jghost;
        int jbeg,jend;
        if(side == left) {
            jghost = data->beg[JDIR];
            jbeg = 0;
            jend = data->beg[JDIR];
            //return;
        }
        else if(side==right) {
            jghost = data->end[JDIR]-1;
            jbeg=data->end[JDIR];
            jend=data->np_tot[JDIR];
        }


        idefix_for("UserDefBoundary",
          0, data->np_tot[KDIR],
          jbeg, jend,
          0, data->np_tot[IDIR],
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        Vc(RHO,k,j,i) = Vc(RHO,k,jghost,i);
                        Vc(PRS,k,j,i) = Vc(PRS,k,jghost,i);
                        Vc(VX1,k,j,i) = ZERO_F;
                        Vc(VX2,k,j,i) = ZERO_F;
                        Vc(VX3,k,j,i) = ZERO_F;

                    });

        idefix_for("UserDefBoundary_X1s",
          0, data->np_tot[KDIR],
          jbeg, jend,
          0, data->np_tot[IDIR]+1,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        Vs(BX1s,k,j,i) = ZERO_F;
                    });

        idefix_for("UserDefBoundary_X3s",
          0, data->np_tot[KDIR]+1,
          jbeg, jend,
          0, data->np_tot[IDIR],
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        Vs(BX3s,k,j,i) = ZERO_F;
                    });


    }

}

// Default constructor

void ComputeUserVars(DataBlock & data, UserDefVariablesContainer &variables) {

  // Make references to the user-defined arrays (variables is a container of IdefixHostArray3D)
  // Note that the labels should match the variable names in the input file
  IdefixHostArray3D<real> Er  = variables["Er"];
  IdefixHostArray3D<real> Eth = variables["Eth"];

    Kokkos::deep_copy(Er,data.hydro->emf->Ex1);
    Kokkos::deep_copy(Eth,data.hydro->emf->Ex2);


}

// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  // Set the function for userdefboundary
  data.hydro->EnrollUserDefBoundary(&UserdefBoundary);
  data.hydro->EnrollUserSourceTerm(&MySourceTerm);
  data.hydro->EnrollInternalBoundary(&InternalBoundary);
  if(data.haveFargo)
    data.fargo->EnrollVelocity(&FargoVelocity);
  output.EnrollUserDefVariables(&ComputeUserVars);

  gammaGlob=data.hydro->eos->GetGamma();
  epsilonGlob = input.Get<real>("Setup","epsilon",0);
  densityFloorGlob = input.Get<real>("Setup","densityFloor",0);
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

    real epsilon=0.1;
    real beta=1000;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                real r=d.x[IDIR](i);
                real th=d.x[JDIR](j);
                real R=r*sin(th);
                real z=r*cos(th);
                real Vk=1.0/pow(R,0.5);

                real cs2=(epsilon*Vk)*(epsilon*Vk);

                d.Vc(RHO,k,j,i) = 1.0/(R*sqrt(R))*exp(1.0/(cs2)*(1/r-1/R));
                d.Vc(PRS,k,j,i) = cs2*d.Vc(RHO,k,j,i);
                d.Vc(VX1,k,j,i) = 0.0;
                d.Vc(VX2,k,j,i) = 1e-1*(0.5-idfx::randm());
                d.Vc(VX3,k,j,i) = Vk*sqrt(R/r-2.5*cs2);

                d.Vs(BX1s,k,j,i) = 0.0;
                d.Vs(BX2s,k,j,i) = 0.0;
                d.Vs(BX3s,k,j,i) = 0.0;

                real B0 = sqrt(2*cs2/(R*sqrt(R))/sqrt(beta));

                r=d.xl[IDIR](i);
                th=d.xl[JDIR](j);
                R=r*sin(th);
                z=r*cos(th);

                A(IDIR,k,j,i) = 0.0;
                A(JDIR,k,j,i) = 0.0;
                A(KDIR,k,j,i) = B0*epsilon*cos(R/epsilon)*fmax(1-(z*z)/(4*R*R*epsilon*epsilon),ZERO_F);
            }
        }
    }

    // Make the field from the vector potential
    d.MakeVsFromAmag(A);

    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}
