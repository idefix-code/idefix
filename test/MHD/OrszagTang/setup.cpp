#include "idefix.hpp"
#include "setup.hpp"

/*********************************************/
/**
Customized random number generator
Allow one to have consistant random numbers
generators on different architectures.
**/
/*********************************************/

void CheckConservation(DataBlock &data) {
  static real firstCall{true};
  static std::array<real,DefaultPhysics::nvar> consArray;
  data.hydro->ConvertPrimToCons();
  auto Uc = data.hydro->Uc;
  auto dV = data.dV;
  //idfx::cout << "Analysis: checking conserved quantities..." << std::endl;
  #ifdef SINGLE_PRECISION
    const real threshold = 1e-4;
  #else
    const real threshold = 1e-13;
  #endif
  for(int nv = 0 ; nv < DefaultPhysics::nvar ; nv++) {
    real cons = 0;
    idefix_reduce("Conserved quantity reduction",
        data.beg[KDIR], data.end[KDIR],
        data.beg[JDIR], data.end[JDIR],
        data.beg[IDIR], data.end[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i, real &c) {
                c += dV(k,j,i)*Uc(nv,k,j,i);
            },
        Kokkos::Sum<real>(cons));
    #ifdef WITH_MPI
      MPI_Allreduce(MPI_IN_PLACE, &cons, 1, realMPI, MPI_SUM, MPI_COMM_WORLD);
    #endif
    if(firstCall) {
      consArray[nv] = cons;
    } else {
      real err = std::fabs((consArray[nv]-cons));
      if(err>threshold) {
        std::stringstream str;
        str << "Quantity " << data.hydro->VcName[nv] << " is not conserved" << std::endl;
        std::cout << "Error on " << data.hydro->VcName[nv] << " is " << err << std::endl;
        std::cout << "Original=" << consArray[nv] << " New=" << cons << std::endl;
        IDEFIX_ERROR(str);
      }
    }
  }
  firstCall=false;
  //idfx::cout << "Analysis: done." << std::endl;
}
// Default constructor

// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  output.EnrollAnalysis(&CheckConservation);
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);

    bool haveTracer = data.hydro->haveTracer;

    real B0=1.0/sqrt(4.0*M_PI);

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                real x=d.x[IDIR](i);
                real y=d.x[JDIR](j);

                d.Vc(RHO,k,j,i) = 25.0/(36.0*M_PI);
                d.Vc(PRS,k,j,i) = 5.0/(12.0*M_PI);
                d.Vc(VX1,k,j,i) = -sin(2.0*M_PI*y);
                d.Vc(VX2,k,j,i) = sin(2.0*M_PI*x);
                #ifdef EVOLVE_VECTOR_POTENTIAL
                  x=d.xl[IDIR](i);
                  y=d.xl[JDIR](j);
                  d.Ve(AX3e,k,j,i) = B0/(2.0*M_PI)*(
                                      cos(2.0*M_PI*y) + cos(4.0*M_PI*x)/2.0);
                #else
                  d.Vs(BX1s,k,j,i) = -B0*sin(2.0*M_PI*y);
                  d.Vs(BX2s,k,j,i) = B0*sin(4.0*M_PI*x);
                #endif
                if(haveTracer) {
                  d.Vc(TRG,k,j,i) = x>0.5?  1.0:0.0;
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
