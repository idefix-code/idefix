#include "idefix.hpp"
#include "setup.hpp"


void FargoVelocity(DataBlock &data, IdefixArray2D<real> &Vphi) {
  IdefixArray1D<real> x1 = data.x[IDIR];
    IdefixArray1D<real> x2 = data.x[JDIR];

  idefix_for("FargoVphi",0,data.np_tot[JDIR], 0, data.np_tot[IDIR],
      KOKKOS_LAMBDA (int j, int i) {
      Vphi(j,i) = 1.0/sqrt(x1(i)*sin(x2(j)));
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

                          Vs(BX2s,k,j,i) = Vs(BX2s,k,j,ighost);


                      });

        idefix_for("UserDefBoundaryX3S",
          0, data->np_tot[KDIR]+1,
          0, data->np_tot[JDIR],
          0, ighost,
                      KOKKOS_LAMBDA (int k, int j, int i) {
                          Vs(BX3s,k,j,i) = Vs(BX3s,k,j,ighost);

            });
    }
}


void Potential(DataBlock& data, const real t, IdefixArray1D<real>& x1, IdefixArray1D<real>& x2, IdefixArray1D<real>& x3, IdefixArray3D<real>& phi) {
    // To simplify, we solve here for a cylindrical potential, to avoid having to treat the stratification in theta
    idefix_for("Potential",0,data.np_tot[KDIR], 0, data.np_tot[JDIR], 0, data.np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
        phi(k,j,i) = -1.0/(x1(i)*sin(x2(j)));
    });

}



// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  // Set the function for userdefboundary
  data.hydro->EnrollUserDefBoundary(&UserdefBoundary);
  data.gravity->EnrollPotential(&Potential);
  if(data.haveFargo)
    data.fargo->EnrollVelocity(&FargoVelocity);
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);

    // Make vector potential
    #ifndef EVOLVE_VECTOR_POTENTIAL
      IdefixHostArray4D<real> A = IdefixHostArray4D<real>("Setup_VectorPotential", 3, data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);
    #endif

    real B0=1.0e-6;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                real r=d.x[IDIR](i);
                real th=d.x[JDIR](j);
                real R=r*sin(th);
                real z=r*cos(th);
                real Vk=1.0/pow(R,0.5);

                d.Vc(RHO,k,j,i) = 1.0;
                d.Vc(VX1,k,j,i) = 0.0;
                d.Vc(VX2,k,j,i) = 0.0;
                d.Vc(VX3,k,j,i) = Vk;

                d.Vs(BX1s,k,j,i) = 0.0;
                d.Vs(BX2s,k,j,i) = 0.0;
                d.Vs(BX3s,k,j,i) = 0.0;

                r=d.xl[IDIR](i);
                th=d.xl[JDIR](j);
                R=r*sin(th);
                z=r*cos(th);
                #ifdef EVOLVE_VECTOR_POTENTIAL
                  d.Ve(AX1e,k,j,i) = 0.0;
                  d.Ve(AX2e,k,j,i) = 0.0;
                  d.Ve(AX3e,k,j,i) = B0*cos(6*d.x[KDIR](k));
                #else
                  A(IDIR,k,j,i) = 0.0;
                  A(JDIR,k,j,i) = 0.0;
                  A(KDIR,k,j,i) = B0*cos(6*d.x[KDIR](k));
                #endif
            }
        }
    }

    #ifndef EVOLVE_VECTOR_POTENTIAL
      // Make the field from the vector potential
      d.MakeVsFromAmag(A);
    #endif

    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}
