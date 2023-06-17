#include "idefix.hpp"
#include "setup.hpp"


real epsilonGlob;
real betaGlob;
real HidealGlob;
real AmMidGlob;
real gammaGlob;
real densityFloorGlob;


void MySoundSpeed(DataBlock &data, const real t, IdefixArray3D<real> &cs) {
  IdefixArray1D<real> r=data.x[IDIR];
  IdefixArray1D<real> th=data.x[JDIR];
  real epsilon = epsilonGlob;
  idefix_for("MySoundSpeed",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real R = r(i)*sin(th(j));
                cs(k,j,i) = epsilon/sqrt(R);
              });
}


// User-defined boundaries
void UserdefBoundary(Hydro *hydro, int dir, BoundarySide side, real t) {
    auto *data = hydro->data;

    if( (dir==IDIR) && (side == left)) {
        IdefixArray4D<real> Vc = hydro->Vc;
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
                        if(Vc(VX1,k,j,ighost)>=ZERO_F) Vc(VX1,k,j,i) = - Vc(VX1,k,j,2*ighost-i);
                               else Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost);
                        Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost);
                        Vc(VX3,k,j,i) = R*Omega;
                    });
    }

    if( dir==JDIR) {
        IdefixArray4D<real> Vc = hydro->Vc;
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
                        Vc(VX1,k,j,i) = ZERO_F;
                        Vc(VX2,k,j,i) = ZERO_F;
                        Vc(VX3,k,j,i) = ZERO_F;
                    });
    }

}


// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  // Set the function for userdefboundary
  data.hydro->EnrollUserDefBoundary(&UserdefBoundary);
  data.hydro->EnrollIsoSoundSpeed(&MySoundSpeed);
  epsilonGlob = input.Get<real>("Setup","epsilon",0);
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);
    real r,th;
    real epsilon=epsilonGlob;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                r=d.x[IDIR](i);
                th=d.x[JDIR](j);
                real R=r*sin(th);
                real z=r*cos(th);
                real Vk=1.0/pow(R,0.5);

                real cs2=(epsilon*Vk)*(epsilon*Vk);

                d.Vc(RHO,k,j,i) = 1.0/(R*sqrt(R))*exp(1.0/(cs2)*(1/r-1/R));
                d.Vc(VX1,k,j,i) = 0.0;
                d.Vc(VX2,k,j,i) = 1e-1*(0.5-idfx::randm());
                d.Vc(VX3,k,j,i) = Vk*sqrt(R/r-2.5*epsilon*epsilon);

            }
        }
    }


    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}
