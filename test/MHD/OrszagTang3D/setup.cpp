#include "idefix.hpp"
#include "setup.hpp"

/*********************************************/
/**
Customized random number generator
Allow one to have consistant random numbers
generators on different architectures.
**/
/*********************************************/

Output* myOutput;
int outnum;

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
    real x,y,z;
    IdefixHostArray4D<real> Ve;

    #ifndef EVOLVE_VECTOR_POTENTIAL
    Ve = IdefixHostArray4D<real>("Potential vector",3, d.np_tot[KDIR]+1, d.np_tot[JDIR]+1, d.np_tot[IDIR]+1);
    #else
    Ve = d.Ve;
    #endif

    bool haveTracer = data.hydro->haveTracer;

    real B0=1.0/sqrt(4.0*M_PI);

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                x=d.x[IDIR](i);
                y=d.x[JDIR](j);
                z=d.x[KDIR](k);

                d.Vc(RHO,k,j,i) = 25.0/(36.0*M_PI);
                d.Vc(PRS,k,j,i) = 5.0/(12.0*M_PI);
                d.Vc(VX1,k,j,i) = -sin(2.0*M_PI*y);
                d.Vc(VX2,k,j,i) = sin(2.0*M_PI*x)+cos(2.0*M_PI*z);
                d.Vc(VX3,k,j,i) = cos(2.0*M_PI*x);

                real xl=d.xl[IDIR](i);
                real yl=d.xl[JDIR](j);
                real zl=d.xl[KDIR](k);
                Ve(IDIR,k,j,i) = B0/(2.0*M_PI)*(cos(2.0*M_PI*yl));
                Ve(JDIR,k,j,i) = B0/(2.0*M_PI)*sin(2.0*M_PI*xl);
                Ve(KDIR,k,j,i) = B0/(2.0*M_PI)*(
                                    cos(2.0*M_PI*yl) + cos(4.0*M_PI*xl)/2.0);

                if(haveTracer) {
                  d.Vc(TRG  ,k,j,i) = x>0.5?  1.0:0.0;
                  d.Vc(TRG+1,k,j,i) = z>0.5?  1.0:0.0;
                }
            }
        }
    }

    #ifndef EVOLVE_VECTOR_POTENTIAL
    d.MakeVsFromAmag(Ve);
    #endif
    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}



// Do a specifically designed user step in the middle of the integration
void ComputeUserStep(DataBlock &data, real t, real dt) {

}
