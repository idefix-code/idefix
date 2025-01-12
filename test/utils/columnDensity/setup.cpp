#include "idefix.hpp"
#include "setup.hpp"
#include "column.hpp"


Column *columnX1Left;
Column *columnX1Right;
Column *columnX2Left;
Column *columnX2Right;
Column *columnX3Left;
Column *columnX3Right;

// Analyse data to check that column density works as expected
void Analysis(DataBlock & data) {

  DataBlockHost d(data);

  columnX1Left->ComputeColumn(data.hydro->Vc);
  columnX1Right->ComputeColumn(data.hydro->Vc);
  columnX2Left->ComputeColumn(data.hydro->Vc);
  columnX2Right->ComputeColumn(data.hydro->Vc);
  columnX3Left->ComputeColumn(data.hydro->Vc);
  columnX3Right->ComputeColumn(data.hydro->Vc);

  IdefixArray3D<real> columnDensityLeft, columnDensityRight;
  IdefixArray3D<real>::HostMirror columnDensityLeftHost, columnDensityRightHost;
  // IDIR
  columnDensityLeft = columnX1Left->GetColumn();
  columnDensityRight = columnX1Right->GetColumn();
  columnDensityLeftHost = Kokkos::create_mirror_view(columnDensityLeft);
  columnDensityRightHost = Kokkos::create_mirror_view(columnDensityRight);
  Kokkos::deep_copy(columnDensityLeftHost,columnDensityLeft);
  Kokkos::deep_copy(columnDensityRightHost,columnDensityRight);

  real errMax = 0.0;

  // Because this particular setup has a constant density =1, we expect
  // the column density from the left to match the cell right coordinates,
  // and the column density from the right to match L-the cell left coordinate (where L is the box size)
  // here we have L=1
  // This routine checks that all of the available column densities do match the expected values.
  for(int k = data.beg[KDIR]; k < data.end[KDIR] ; k++) {
    for(int j = data.beg[JDIR]; j < data.end[JDIR] ; j++) {
      for(int i = data.beg[IDIR]; i < data.end[IDIR] ; i++) {
        real err = std::fabs(columnDensityLeftHost(k,j,i)-d.xr[IDIR](i));
        if(err>errMax) errMax=err;
        err = std::fabs(columnDensityRightHost(k,j,i)-(1-d.xl[IDIR](i)));
        if(err>errMax) errMax=err;
      }
    }
  }
  idfx::cout << "Error on column density in IDIR=" << std::scientific << errMax << std::endl;
  if(errMax>1e-14) {
    IDEFIX_ERROR("Error above tolerance");
  }
  // JDIR
  columnDensityLeft = columnX2Left->GetColumn();
  columnDensityRight = columnX2Right->GetColumn();
  columnDensityLeftHost = Kokkos::create_mirror_view(columnDensityLeft);
  columnDensityRightHost = Kokkos::create_mirror_view(columnDensityRight);
  Kokkos::deep_copy(columnDensityLeftHost,columnDensityLeft);
  Kokkos::deep_copy(columnDensityRightHost,columnDensityRight);

  errMax = 0.0;

  for(int k = data.beg[KDIR]; k < data.end[KDIR] ; k++) {
    for(int j = data.beg[JDIR]; j < data.end[JDIR] ; j++) {
      for(int i = data.beg[IDIR]; i < data.end[IDIR] ; i++) {
        real err = std::fabs(columnDensityLeftHost(k,j,i)-d.xr[JDIR](j));
        if(err>errMax) errMax=err;
        err = std::fabs(columnDensityRightHost(k,j,i)-(1-d.xl[JDIR](j)));
        if(err>errMax) errMax=err;
      }
    }
  }
  idfx::cout << "Error on column density in JDIR=" << std::scientific << errMax << std::endl;
  if(errMax>1e-14) {
    IDEFIX_ERROR("Error above tolerance");
  }
  // KDIR
  columnDensityLeft = columnX3Left->GetColumn();
  columnDensityRight = columnX3Right->GetColumn();
  columnDensityLeftHost = Kokkos::create_mirror_view(columnDensityLeft);
  columnDensityRightHost = Kokkos::create_mirror_view(columnDensityRight);
  Kokkos::deep_copy(columnDensityLeftHost,columnDensityLeft);
  Kokkos::deep_copy(columnDensityRightHost,columnDensityRight);

  errMax = 0.0;

  for(int k = data.beg[KDIR]; k < data.end[KDIR] ; k++) {
    for(int j = data.beg[JDIR]; j < data.end[JDIR] ; j++) {
      for(int i = data.beg[IDIR]; i < data.end[IDIR] ; i++) {
        real err = std::fabs(columnDensityLeftHost(k,j,i)-d.xr[KDIR](k));
        if(err>errMax) errMax=err;
        err = std::fabs(columnDensityRightHost(k,j,i)-(1-d.xl[KDIR](k)));
        if(err>errMax) errMax=err;
      }
    }
  }
  idfx::cout << "Error on column density in KDIR=" << std::scientific << errMax << std::endl;
  if(errMax>1e-14) {
    IDEFIX_ERROR("Error above tolerance");
  }

}

void InternalBoundary(Fluid<DefaultPhysics> * hydro, const real t) {
  IdefixArray4D<real> Vc = hydro->Vc;
  idefix_for("InternalBoundary",0,hydro->data->np_tot[KDIR],
                                0,hydro->data->np_tot[JDIR],
                                0,hydro->data->np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                // Cancel any motion that could be happening
                Vc(VX1,k,j,i) = 0.0;
                Vc(VX2,k,j,i) = 0.0;
                Vc(VX3,k,j,i) = 0.0;
              });
}


// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  output.EnrollAnalysis(&Analysis);
  data.hydro->EnrollInternalBoundary(&InternalBoundary);

  columnX1Left = new Column(IDIR, 1, RHO, &data);
  columnX1Right = new Column(IDIR, -1, RHO, &data);
  columnX2Left = new Column(JDIR, 1, RHO, &data);
  columnX2Right = new Column(JDIR, -1, RHO, &data);
  columnX3Left = new Column(KDIR, 1, RHO, &data);
  columnX3Right = new Column(KDIR, -1, RHO, &data);
  // Initialise the output file
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

                d.Vc(RHO,k,j,i) = 1.0;
                d.Vc(VX1,k,j,i) = ZERO_F;
                d.Vc(PRS,k,j,i) = 1.0;

            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}
