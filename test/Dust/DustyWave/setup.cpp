#include "idefix.hpp"
#include "setup.hpp"

#define  FILENAME    "timevol.dat"

// Analyse data to produce an output
void Analysis(DataBlock & data) {

  auto Vc=data.hydro->Vc;
  real Ek, Erho;
  idefix_reduce("Vx2",data.beg[KDIR],data.end[KDIR],data.beg[JDIR],data.end[JDIR],data.beg[IDIR],data.end[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i, real &vx2) {
                vx2 += Vc(VX1,k,j,i) * Vc(VX1,k,j,i);
              }, Kokkos::Sum<real>(Ek) );

  idefix_reduce("RHO2",data.beg[KDIR],data.end[KDIR],data.beg[JDIR],data.end[JDIR],data.beg[IDIR],data.end[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i, real &rho2) {
                rho2 += (Vc(RHO,k,j,i)-1)* (Vc(RHO,k,j,i)-1);
              }, Kokkos::Sum<real>(Erho) );

  Ek /= (data.np_int[KDIR]*data.np_int[JDIR]*data.np_int[IDIR]);
  Erho /= (data.np_int[KDIR]*data.np_int[JDIR]*data.np_int[IDIR]);

  if(idfx::prank == 0) {
    std::ofstream f;
    f.open(FILENAME,std::ios::app);
    f.precision(10);
    f << std::scientific << data.t << "\t" << Ek << "\t" << Erho << std::endl;
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
      f << "t\t\t vx2 \t\t rho2" << std::endl;
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


    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {

                d.Vc(RHO,k,j,i) = 1.0;
                d.dustVc[0](RHO,k,j,i) = 1.0;

                d.Vc(VX1,k,j,i) = 0.01*sin(2.0*M_PI*d.x[IDIR](i));
                d.dustVc[0](VX1,k,j,i) = 0.0;

#if HAVE_ENERGY
                d.Vc(PRS,k,j,i) = 1.0;
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
