// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// This setup is the base for a wind-driven accretion disc dominated by Ambipolar diffusion


#include "idefix.hpp"
#include "setup.hpp"

real RinGlob;
real sigma0Glob;
real sigmaSlopeGlob;
real h0Glob;

real h0TopGlob;
real dampingZoneGlob;
real tauDampGlob;
real betaGlob;
real HidealGlob;
real AmMidGlob;
real gammaGlob;
real densityFloorGlob;


void Ambipolar(DataBlock& data, real t, IdefixArray3D<real> &xAin) {
  IdefixArray3D<real> xA = xAin;
  IdefixArray1D<real> x1=data.x[IDIR];
  IdefixArray1D<real> x3=data.x[KDIR];
  IdefixArray4D<real> Vc=data.hydro.Vc;

  real Hideal = HidealGlob;
  real h0 = h0Glob;
  real AmMid = AmMidGlob;

  idefix_for("Ambipolar",0,data.np_tot[KDIR],
                         0,data.np_tot[JDIR],
                         0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real z=x3(k);
                real R=x1(i);
                real Omega=pow(R,-1.5);

                real q = z/(Hideal*R*h0);
                real Am = AmMid * exp(q*q*q*q);

                real etamax = 10*Omega *h0*h0*R*R;  // Rm=10 cap

                real B2 = Vc(BX1,k,j,i)*Vc(BX1,k,j,i)
                          +Vc(BX2,k,j,i)*Vc(BX2,k,j,i)
                          +Vc(BX3,k,j,i)*Vc(BX3,k,j,i);
                real eta = B2/(Omega*Am*Vc(RHO,k,j,i));
                if(eta>etamax)
                  xA(k,j,i) = etamax/B2;
                else
                  xA(k,j,i) = 1.0/(Omega*Am*Vc(RHO,k,j,i));
              });
}

void Resistivity(DataBlock& data, real t, IdefixArray3D<real> &etain) {
  IdefixArray3D<real> eta = etain;
  IdefixArray1D<real> x1=data.x[IDIR];
  IdefixArray4D<real> Vc=data.hydro.Vc;

  real h0 = h0Glob;
  real Rin = RinGlob;
  real etaMax = h0*h0*Rin*Rin*pow(Rin,-1.5);
  real dampingZone=dampingZoneGlob;

  idefix_for("Resistivity",0,data.np_tot[KDIR],
                           0,data.np_tot[JDIR],
                           0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real R=x1(i);

                if(R<dampingZone*Rin) {
                  eta(k,j,i) = etaMax;
                } else {
                  eta(k,j,i) = ZERO_F;
                }
              });
}

//loi de temperature dans le disque transition entre midplane et corona
// (en tan hyperbolique) car plus chaud dans la couronne -> rajouter ce chauffage
void MySourceTerm(DataBlock &data, const real t, const real dtin) {
  IdefixArray4D<real> Vc = data.hydro.Vc;
  IdefixArray4D<real> Uc = data.hydro.Uc;
  IdefixArray1D<real> x1=data.x[IDIR];
  IdefixArray1D<real> x3=data.x[KDIR];
  real Rin = RinGlob;
  real h0Top = h0TopGlob;
  real h0 = h0Glob;
  real dampingZone=dampingZoneGlob;
  real tauDamp=tauDampGlob;
  real tauGlob=0.1;
  real gamma_m1=gammaGlob-1.0;
  real dt=dtin;
  real Hideal=HidealGlob;
  idefix_for("MySourceTerm",0,data.np_tot[KDIR],
                            0,data.np_tot[JDIR],
                            0,data.np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      real R=x1(i);
      real z=x3(k);
      real cs2, tau, taud;
      real Zh = FABS(z/R)/h0;
      real csdisk = h0/sqrt(R);
      real cscorona = h0Top/sqrt(R);
      cs2=0.5*(csdisk*csdisk+cscorona*cscorona)
         +0.5*(cscorona*cscorona-csdisk*csdisk)*tanh(6*log(Zh/Hideal));
      tau=tauGlob*sqrt(R); // 0.1 temps orbital local (implementation du beta cooling).
                            // Si isotherme, tps cooling = 0, et VSI violente
                            //-> ici VSI amortie.
      taud=tauDamp*sqrt(R); // tauDamp temps orbitaux locaux.

      // Cooling /heatig function
      real Ptarget = cs2*Vc(RHO,k,j,i);

      // temperature relaxe vers un profil de temperature defini par cs2,
      // ENG energie totale
      Uc(ENG,k,j,i) += -dt*(Vc(PRS,k,j,i)-Ptarget)/(tau*gamma_m1);

    // Velocity relaxation -> amortir les mouvements pour R<dampingZone*Rin
    // (like WKZ), MX1 = RHO*VX1
    if(R<dampingZone*Rin) {
        Uc(MX1,k,j,i) += -dt*(Vc(VX1,k,j,i)*Vc(RHO,k,j,i))/taud;
        Uc(MX3,k,j,i) += -dt*(Vc(VX3,k,j,i)*Vc(RHO,k,j,i))/taud;
    }
  });
}

// detecte quand la vitesse alfven depasse un seuil,
// modification de la densité localement (limiter la vitesse d'alfven)
// force la densité dans la couronne à ne pas etre trop petite
void InternalBoundary(DataBlock& data, const real t) {
  IdefixArray4D<real> Vc = data.hydro.Vc;
  IdefixArray4D<real> Vs = data.hydro.Vs;

  real vAmax=1.0*pow(RinGlob,-0.5);

  real densityFloor = densityFloorGlob;
  idefix_for("InternalBoundary",0,data.np_tot[KDIR],
                                0,data.np_tot[JDIR],
                                0,data.np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      real b2=EXPAND(Vc(BX1,k,j,i)*Vc(BX1,k,j,i) ,
                    +Vc(BX2,k,j,i)*Vc(BX2,k,j,i) ,
                    +Vc(BX3,k,j,i)*Vc(BX3,k,j,i));

      real va2=b2/Vc(RHO,k,j,i);
      real myMax=vAmax;
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
    });
}

// User-defined boundaries
void UserdefBoundary(DataBlock& data, int dir, BoundarySide side, real t) {
  // direction R (min)
  if( (dir==IDIR) && (side == left)) {
    IdefixArray4D<real> Vc = data.hydro.Vc;
    IdefixArray4D<real> Vs = data.hydro.Vs;
    IdefixArray1D<real> x1 = data.x[IDIR];
    IdefixArray1D<real> x3 = data.x[KDIR];
    real Rin = RinGlob;
      real h0 = h0Glob;
    real h0Top = h0TopGlob;
    int ighost = data.nghost[IDIR];
    // real Omegain=1.0; // if solid rotation
    real Omegain=pow(Rin,-1.5); // if keplerian rotation
    real Hideal=HidealGlob;
    idefix_for("UserDefBoundary",0,data.np_tot[KDIR],
                                 0,data.np_tot[JDIR],
                                 0,ighost,
      KOKKOS_LAMBDA (int k, int j, int i) {
        real R=x1(i);
        real z=x3(k);

        real Zh = FABS(z/R)/h0;
        real csdisk = h0/sqrt(R);
        real cscorona = h0Top/sqrt(R);
        real cs2=0.5*(csdisk*csdisk+cscorona*cscorona)
                +0.5*(cscorona*cscorona-csdisk*csdisk)*tanh(6*log(Zh/Hideal));

        Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost);
        Vc(PRS,k,j,i) = cs2*Vc(RHO,k,j,i);

        Vc(VX1,k,j,i) = - Vc(VX1,k,j,2*ighost-i-1); // les flux sont ~ nuls a l'interface
        Vc(VX2,k,j,i) = pow(R,-0.5); // rotation solide -> si keplerien : pow(r,-0.5)
        Vc(VX3,k,j,i) = ZERO_F; // vitesse verticale nulle
    });

    idefix_for("UserDefBoundary",0,data.np_tot[KDIR],
                                 0,data.np_tot[JDIR]+JOFFSET,
                                 0,ighost,
      KOKKOS_LAMBDA (int k, int j, int i) {
          Vs(BX2s,k,j,i) = Vs(BX2s,k,j,ighost) ;
    });

    idefix_for("UserDefBoundary",0,data.np_tot[KDIR]+KOFFSET,
                                 0,data.np_tot[JDIR],
                                 0,ighost,
      KOKKOS_LAMBDA (int k, int j, int i) {
          //Vs(BX2s,k,j,i) = Vs(BX2s,k,j,ighost)*pow(x1(i)/x1(ighost),-5.0/4.0) ;
          Vs(BX3s,k,j,i) = Vs(BX3s,k,j,ighost);
    });
  }

    // Self-similar outside

    if( (dir==IDIR) && (side == right)) {
      IdefixArray4D<real> Vc = data.hydro.Vc;
      IdefixArray4D<real> Vs = data.hydro.Vs;
      IdefixArray1D<real> x1 = data.x[IDIR];

      int ighost = data.end[IDIR]-1;

      idefix_for("UserDefBoundary",0,data.np_tot[KDIR],
                                   0,data.np_tot[JDIR],
                                   data.end[IDIR],data.np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
            Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost)*pow(x1(i)/x1(ighost),-1.5);
            Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost)*pow(x1(i)/x1(ighost),-2.5);
            Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost)*pow(x1(i)/x1(ighost),-0.5);;
            Vc(VX3,k,j,i) = Vc(VX3,k,j,ighost)*pow(x1(i)/x1(ighost),-0.5);;
            Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost)*pow(x1(i)/x1(ighost),-0.5);;
        });
      idefix_for("UserDefBoundary",0,data.np_tot[KDIR],
                                   0,data.np_tot[JDIR]+JOFFSET,
                                   data.end[IDIR],data.np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vs(BX2s,k,j,i) = Vs(BX2s,k,j,ighost)*pow(x1(i)/x1(ighost),-5.0/4.0);
        });

      idefix_for("UserDefBoundary",0,data.np_tot[KDIR]+KOFFSET,
                                   0,data.np_tot[JDIR],
                                   data.end[IDIR],data.np_tot[IDIR],
      KOKKOS_LAMBDA (int k, int j, int i) {
          Vs(BX3s,k,j,i) = Vs(BX3s,k,j,ighost)*pow(x1(i)/x1(ighost),-5.0/4.0) ;
      });
    }
}

void Potential(DataBlock& data, const real t, IdefixArray1D<real>& x1,
                                              IdefixArray1D<real>& x2,
                                              IdefixArray1D<real>& x3,
                                              IdefixArray3D<real>& phi) {
    idefix_for("Potential",0,data.np_tot[KDIR],
                           0, data.np_tot[JDIR],
                           0, data.np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
        phi(k,j,i) = -1.0/sqrt(x1(i)*x1(i)+x3(k)*x3(k));
    });
}


// Default constructor

// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
    // Set the function for userdefboundary
    data.hydro.EnrollUserDefBoundary(&UserdefBoundary);
    data.hydro.EnrollGravPotential(&Potential);
    data.hydro.EnrollAmbipolarDiffusivity(&Ambipolar);
    data.hydro.EnrollOhmicDiffusivity(&Resistivity);
    data.hydro.EnrollUserSourceTerm(&MySourceTerm);
    data.hydro.EnrollInternalBoundary(&InternalBoundary);
    gammaGlob=data.hydro.GetGamma();
    RinGlob = grid.xbeg[IDIR];
    sigma0Glob = input.GetReal("Setup","sigma0",0);
    sigmaSlopeGlob = input.GetReal("Setup","sigmaSlope",0);
    h0Glob = input.GetReal("Setup","h0",0);
    // flaringIndexGlob = input.GetReal("Setup","flaringIndex",0);
    h0TopGlob = input.GetReal("Setup","h0Top",0);
    dampingZoneGlob = input.GetReal("Setup","dampingZone",0);
    tauDampGlob = input.GetReal("Setup","tauDamp",0);
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
  IdefixHostArray4D<real> A = IdefixHostArray4D<real>("Setup_VectorPotential",
                                                      3,
                                                      data.np_tot[KDIR],
                                                      data.np_tot[JDIR],
                                                      data.np_tot[IDIR]);
  real Rin=RinGlob;
  real h0=h0Glob;
  real sigma0=sigma0Glob;
  real sigmaSlope=sigmaSlopeGlob;
  real m=-5.0/4.0;
  real B0 = h0*sqrt(2.0/betaGlob);

  for(int k = 0; k < d.np_tot[KDIR] ; k++) {
    for(int j = 0; j < d.np_tot[JDIR] ; j++) {
      for(int i = 0; i < d.np_tot[IDIR] ; i++) {
        real R=d.x[IDIR](i);
        real z=d.x[KDIR](k);
        real Zh = FABS(z/R)/h0;
        real csdisk = h0/sqrt(R);
        real cscorona = h0TopGlob/sqrt(R);
        real cs2=0.5*(csdisk*csdisk+cscorona*cscorona)+
               0.5*(cscorona*cscorona-csdisk*csdisk)*tanh(6*log(Zh/HidealGlob));

        d.Vc(RHO,k,j,i) = sigma0/sqrt(2.0*M_PI)/(h0*R)*pow(R,-sigmaSlope)
                        * exp(1.0/ (csdisk*csdisk)* (1.0/sqrt(R*R+z*z)-1.0/R));
        d.Vc(VX2,k,j,i) = 1.0/sqrt(R)
                        * sqrt( R / sqrt(R*R + z*z) -(1+sigmaSlope)*h0*h0);
        d.Vc(PRS,k,j,i) = cs2*d.Vc(RHO,k,j,i);

        d.Vc(VX1,k,j,i) = ZERO_F;
        d.Vc(VX3,k,j,i) = ZERO_F;

        if(d.Vc(RHO,k,j,i) < densityFloorGlob) {
          real T2=d.Vc(PRS,k,j,i)/d.Vc(RHO,k,j,i);
          d.Vc(RHO,k,j,i) = densityFloorGlob;
          d.Vc(PRS,k,j,i) = T2*d.Vc(RHO,k,j,i);
        }

        // Vector potential on the corner
        R=d.xl[IDIR](i);

        A(IDIR,k,j,i) = ZERO_F;
        A(KDIR,k,j,i) = ZERO_F;
        A(JDIR,k,j,i) = B0*(pow(Rin,m+2.0)/R * (-1.0/(m+2.0))
                          + pow(R,m+1.0)/(m+2.0));
      }
    }
  }

  // Make the field from the vector potential
  d.MakeVsFromAmag(A);


  // Send it all, if needed
  d.SyncToDevice();
}






// Do a specifically designed user step in the middle of the integration
void ComputeUserStep(DataBlock &data, real t, real dt) {

}
