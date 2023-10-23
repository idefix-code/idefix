#include "idefix.hpp"
#include "setup.hpp"


#define FILENAME  "analysis.dat"
real amplitude;
// Analyse data to get the amplitude of the k=1 temperature perturbation
void Analysis(DataBlock & data) {
  DataBlockHost d(data);
  d.SyncFromDevice();
  real Tmode = 0;
  for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
    real T = d.Vc(PRS,d.beg[KDIR],d.beg[JDIR],i) / d.Vc(RHO,d.beg[KDIR],d.beg[JDIR],i);
    Tmode += T * sin(2.0*M_PI*d.x[IDIR](i));
  }
  Tmode = 2*Tmode/d.np_int[JDIR];

  std::ofstream f;
  f.open(FILENAME,std::ios::app);
  f.precision(10);
  f << std::scientific << data.t << "\t" << Tmode << std::endl;
  f.close();

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

  amplitude = input.Get<real>("Setup","amplitude",0);
  // Initialise the output file
  std::ofstream f;
  f.open(FILENAME,std::ios::trunc);
  f << "t\t\t T" << std::endl;
  f.close();
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

                d.Vc(RHO,k,j,i) = 1 - amplitude*sin(2.0*M_PI*d.x[IDIR](i));
                d.Vc(VX1,k,j,i) = ZERO_F;
                d.Vc(PRS,k,j,i) = 1.0;

            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}
