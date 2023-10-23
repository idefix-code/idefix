#include "idefix.hpp"
#include "setup.hpp"

/*********************************************/
/**
Customized random number generator
Allow one to have consistant random numbers
generators on different architectures.
**/
/*********************************************/

// User-defined boundaries
void UserdefBoundary(Hydro *hydro, int dir, BoundarySide side, real t) {
    auto *data = hydro->data;
    if( (dir==IDIR) && (side == left)) {
        IdefixArray4D<real> Vc = hydro->Vc;
        IdefixArray4D<real> Vs = hydro->Vs;
        IdefixArray1D<real> x1 = data->x[IDIR];

        int ighost = data->nghost[IDIR];
        idefix_for("UserDefBoundary",
          0, data->np_tot[KDIR],
          0, data->np_tot[JDIR],
          0, ighost,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost);
                        Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost);
                        Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost) * x1(i)/x1(ighost);
                        Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost) * x1(i)/x1(ighost);

                        Vs(BX2s,k,j,i) = Vs(BX2s,k,j,ighost);

                    });
    }

}

// Analyse data to produce an output
void ComputeUserVars(DataBlock & data, UserDefVariablesContainer &variables) {
  DataBlockHost dHost(data);

  IdefixHostArray3D<real> vx = variables["Vx"];
  IdefixHostArray3D<real> vy = variables["Vy"];

  for(int k = 0; k < dHost.np_tot[KDIR] ; k++) {
    for(int j = 0; j < dHost.np_tot[JDIR] ; j++) {
      for(int i = 0; i < dHost.np_tot[IDIR] ; i++) {
        real r = dHost.x[IDIR](i);
        real theta = dHost.x[JDIR](j);
        real x = r*cos(theta);
        real y = r*sin(theta);

        vx(k,j,i) = (dHost.Vc(VX1,k,j,i)*x - dHost.Vc(VX2,k,j,i)*y)/r;
        vy(k,j,i) = (dHost.Vc(VX1,k,j,i)*y + dHost.Vc(VX2,k,j,i)*x)/r;
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
    output.EnrollUserDefVariables(&ComputeUserVars);

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

    real r0 = 0.1;
    real r1 = 0.115;
    real omega = 20;
    real B0 = 5.0/sqrt(4.0*M_PI);
    real vphi,f,r,th;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                r=d.x[IDIR](i);
                th=d.x[JDIR](j);
                real y=r*sin(th);

                f=(r1-r)/(r1-r0);

                if(r<r0) {
                    d.Vc(RHO,k,j,i) = 10.0;
                    vphi = omega*r;
                }
                else if(r<=r1) {
                    d.Vc(RHO,k,j,i) = 1.0+9.0*f;
                    vphi = f*omega*r0;
                }
                else {
                    d.Vc(RHO,k,j,i) = 1.0;
                    vphi = 0.0;
                }
                d.Vc(PRS,k,j,i) = ONE_F;
                d.Vc(VX1,k,j,i) = 0.0;
                d.Vc(VX2,k,j,i) = vphi;



                A(IDIR,k,j,i) = 0.0;
                A(JDIR,k,j,i) = 0.0;
                A(KDIR,k,j,i) = -y*B0;
            }
        }
    }

    // Make the field from the vector potential
    d.MakeVsFromAmag(A);


    // Send it all, if needed
    d.SyncToDevice();
}
