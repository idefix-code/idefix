#include "idefix.hpp"
#include "setup.hpp"

#define  FILENAME    "timevol.dat"

void MyDrag(DataBlock *data, int n, real beta, IdefixArray3D<real> &gamma) {
  // Compute the drag coefficient gamma from the input beta
  auto VcGas = data->hydro->Vc;
  auto VcDust = data->dust[n]->Vc;

  idefix_for("MyDrag",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      gamma(k,j,i) = beta/VcGas(RHO,k,j,i)/VcDust(RHO,k,j,i);
    });
}


void ApplyBoundary(DataBlock *data, IdefixArray4D<real> Vc, int dir, BoundarySide side) {
  if(dir==IDIR) {
    int iref,ibeg,iend;
    if(side == left) {
      iref = data->beg[IDIR];
      ibeg = 0;
      iend = data->beg[IDIR];
    } else {
      iref =data->end[IDIR] - 1;
      ibeg=data->end[IDIR];
      iend=data->np_tot[IDIR];
    }
    int nvar = Vc.extent(0);
    idefix_for("UserDefBoundary",0,data->np_tot[KDIR],0,data->np_tot[JDIR],ibeg,iend,
      KOKKOS_LAMBDA (int k, int j, int i) {
        for(int n = 0 ; n < nvar ; n++) {
          Vc(n,k,j,i) = Vc(n,k,j,iref );
        }
    });
  }
}

void UserdefBoundary(Hydro *hydro, int dir, BoundarySide side, real t) {
  ApplyBoundary(hydro->data, hydro->Vc, dir, side);
}

void UserdefBoundaryDust(Fluid<DustPhysics> *dust, int dir, BoundarySide side, real t) {
  ApplyBoundary(dust->data, dust->Vc, dir, side);
}


// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {

  data.hydro->EnrollUserDefBoundary(&UserdefBoundary);
  if(data.haveDust) {
    int nSpecies = data.dust.size();
    for(int n = 0 ; n < nSpecies ; n++) {
      data.dust[n]->EnrollUserDefBoundary(&UserdefBoundaryDust);
      data.dust[n]->drag->EnrollUserDrag(&MyDrag);
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

    const real x0 = 4.0;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                real x = d.x[IDIR](i);

                d.Vc(RHO,k,j,i) =  (x < x0) ? 1.0 : 16.0;
                d.Vc(VX1,k,j,i) = (x < x0) ? 2.0 : 0.125;

                for(int n = 0 ; n < data.dust.size(); n++) {
                  d.dustVc[n](RHO,k,j,i) = (x < x0) ? 1.0 : 16.0;
                  d.dustVc[n](VX1,k,j,i) = (x < x0) ? 2.0 : 0.125;
                }


            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}
