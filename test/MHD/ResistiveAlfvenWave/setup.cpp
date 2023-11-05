#include "idefix.hpp"
#include "setup.hpp"

#define  FILENAME    "timevol.dat"

// Analyse data to produce an output
void Analysis(DataBlock & data) {
  double etot = 0;
  double etotGlob;
  IdefixArray4D<real> Vc = data.hydro->Vc;

  idefix_reduce("Analysis", data.beg[KDIR],data.end[KDIR],
                data.beg[JDIR],data.end[JDIR],
                data.beg[IDIR],data.end[IDIR],
              KOKKOS_LAMBDA(int k, int j, int i, double &eloc) {
                eloc +=  Vc(VX2,k,j,i)*Vc(VX2,k,j,i)
                        +Vc(VX3,k,j,i)*Vc(VX3,k,j,i)
                        +Vc(BX2,k,j,i)*Vc(BX2,k,j,i)
                        +Vc(BX3,k,j,i)*Vc(BX3,k,j,i);
              }, Kokkos::Sum<double>(etot));

  #ifdef WITH_MPI
    MPI_Reduce(&etot, &etotGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  #else
    etotGlob = etot;
  #endif

  if(idfx::prank == 0) {
    std::ofstream f;
    f.open(FILENAME,std::ios::app);
    f.precision(10);
    f << std::scientific << data.t << "\t" << etotGlob << std::endl;
    f.close();
  }

}


// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {

  output.EnrollAnalysis(&Analysis);

  if(!input.restartRequested) {
    // Initialise the output file
    std::ofstream f;
    f.open(FILENAME,std::ios::trunc);
    f << "t\t\t etot" << std::endl;
    f.close();
  }

}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);

    real B0=1.0;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                real x=d.x[IDIR](i);

                d.Vc(RHO,k,j,i) = 1.0;
                d.Vc(PRS,k,j,i) = 1.0;
                d.Vc(VX1,k,j,i) = 0.0;
#if COMPONENTS >=2
        d.Vc(BX2,k,j,i) = 0.0;
        d.Vc(VX2,k,j,i) = 1.0e-3*sin(2.0*M_PI*x);
#if COMPONENTS == 3
        d.Vc(BX3,k,j,i) = 0.0;
                d.Vc(VX3,k,j,i) = 0.0;
#endif
#endif


                d.Vs(BX1s,k,j,i) = B0;
#if DIMENSIONS >=2
                d.Vs(BX2s,k,j,i) = 0.0;
#if DIMENSIONS >=3
                d.Vs(BX3s,k,j,i) = B0;
#endif
#endif

            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}
