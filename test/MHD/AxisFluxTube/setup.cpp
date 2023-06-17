#include "idefix.hpp"
#include "setup.hpp"



real Rtorus;
real Ztorus;
real Rin;

void ComputeUserVars(DataBlock & data, UserDefVariablesContainer &variables) {
  // Mirror data on Host
  DataBlockHost d(data);

  // Sync it
  d.SyncFromDevice();

  // Make references to the user-defined arrays (variables is a container of IdefixHostArray3D)
  // Note that the labels should match the variable names in the input file
  IdefixHostArray3D<real> divB  = variables["divB"];
  IdefixHostArray3D<real> Er  = variables["Er"];
  IdefixHostArray4D<real> Vs = d.Vs;
  IdefixHostArray3D<real> Ax1 = d.A[IDIR];
  IdefixHostArray3D<real> Ax2 = d.A[JDIR];
  IdefixHostArray3D<real> Ax3 = d.A[KDIR];
  IdefixHostArray3D<real> dV = d.dV;

  for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
    for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
      for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {

        divB(k,j,i) = ((Ax1(k,j,i+1)*Vs(BX1s,k,j,i+1)-Ax1(k,j,i)*Vs(BX1s,k,j,i)) +
                      (Ax2(k,j+1,i)*Vs(BX2s,k,j+1,i)-Ax2(k,j,i)*Vs(BX2s,k,j,i)) +
                      (Ax3(k+1,j,i)*Vs(BX3s,k+1,j,i)-Ax3(k,j,i)*Vs(BX3s,k,j,i)))
                      / dV(k,j,i);

        Er(k,j,i) = d.Ex1(k,j+1,i);
      }
    }
  }
}

void UserdefBoundary(Hydro *hydro, int dir, BoundarySide side, real t) {
    auto *data = hydro->data;
    if( (dir==IDIR) && (side == right)) {
        IdefixArray4D<real> Vc = hydro->Vc;
        IdefixArray4D<real> Vs = hydro->Vs;
        IdefixArray1D<real> x1Arr = data->x[IDIR];
        IdefixArray1D<real> x2Arr = data->x[JDIR];
        IdefixArray1D<real> x3Arr = data->x[KDIR];

        const real sq2 = sqrt(2);

        hydro->boundary->BoundaryFor("UserDefX1",dir,side,
            KOKKOS_LAMBDA (int k, int j, int i) {
                real x1 = x1Arr(i);
                real x2 = x2Arr(j);
                real x3 = x3Arr(k);

                // Vector components in cartesian coordinates

                real ex_r=cos(x3)*sin(x2);
                real ex_t=cos(x3)*cos(x2);
                real ex_p=-sin(x3);

                real ey_r=sin(x3)*sin(x2);
                real ey_t=sin(x3)*cos(x2);
                real ey_p=cos(x3);


                Vc(RHO,k,j,i) = 1.0;
                Vc(PRS,k,j,i) = 1.0;
                Vc(VX1,k,j,i) = (ex_r+ey_r)/sq2;
                Vc(VX2,k,j,i) = (ex_t+ey_t)/sq2;
                Vc(VX3,k,j,i) = (ex_p+ey_p)/sq2;

            });
      hydro->boundary->BoundaryForX2s("UserDefX2s",dir,side,
        KOKKOS_LAMBDA (int k, int j, int i) {
            Vs(BX2s,k,j,i) = ZERO_F;
          });
      hydro->boundary->BoundaryForX3s("UserDefX2s",dir,side,
        KOKKOS_LAMBDA (int k, int j, int i) {
            Vs(BX3s,k,j,i) = ZERO_F;
          });
    }
  }

void Analysis(DataBlock & data) {
  // Mirror data on Host
  data.hydro->boundary->SetBoundaries(data.t);
  data.DumpToFile("analysis");
}

void CoarsenFunction(DataBlock &data) {
  IdefixArray2D<int> coarseningLevel = data.coarseningLevel[KDIR];
  IdefixArray1D<real> th = data.x[JDIR];
  idefix_for("set_coarsening", 0, data.np_tot[JDIR], 0, data.np_tot[IDIR],
      KOKKOS_LAMBDA(int j,int i) {
        coarseningLevel(j,i) = 1;
        if(th(j) < 0.3) {
          coarseningLevel(j,i) = 2;
          /*
          if(th(j)<0.1) {
            coarseningLevel(j,i) = 3;
          }*/
        }
      });
}

// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  output.EnrollUserDefVariables(&ComputeUserVars);
  data.hydro->EnrollUserDefBoundary(&UserdefBoundary);
  output.EnrollAnalysis(&Analysis);
  Rtorus = input.Get<real>("Setup","Rtorus",0);
  Ztorus = input.Get<real>("Setup","Ztorus",0);
  Rin = input.Get<real>("Setup","Rin",0);

  if(data.haveGridCoarsening) {
    data.EnrollGridCoarseningLevels(&CoarsenFunction);
  }
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                real x1=d.x[IDIR](i);
                real x2=d.x[JDIR](j);
                real x3=d.x[KDIR](k);

                real R = x1*sin(x2);      /* Cylindrical radius */
                real Z = x1*cos(x2);      /* Vertical cylindrical coordinate */

                // Vector components in cartesian coordinates

                real ex_r=cos(x3)*sin(x2);
                real ex_t=cos(x3)*cos(x2);
                real ex_p=-sin(x3);

                real ey_r=sin(x3)*sin(x2);
                real ey_t=sin(x3)*cos(x2);
                real ey_p=cos(x3);


                d.Vc(RHO,k,j,i) = 1.0;
                d.Vc(PRS,k,j,i) = 1.0;
                d.Vc(VX1,k,j,i) = (ex_r+ey_r)/sqrt(2);
                d.Vc(VX2,k,j,i) = (ex_t+ey_t)/sqrt(2);
                d.Vc(VX3,k,j,i) = (ex_p+ey_p)/sqrt(2);

                real bphi = 1.0 - (pow(R-Rtorus,2.0) + pow(fabs(Z)-Ztorus,2.0)) / Rin;
                if(bphi<0.0) bphi = 0.0;
                d.Vs(BX1s,k,j,i) = 0.0;
                d.Vs(BX2s,k,j,i) = 0.0;
                d.Vs(BX3s,k,j,i) = 1.0e-20*bphi;


            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}
