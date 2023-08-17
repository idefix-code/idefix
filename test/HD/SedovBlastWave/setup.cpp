#include "idefix.hpp"
#include "setup.hpp"


// Default constructor

real Rstart;
// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  Rstart = input.Get<real>("Setup","Rstart",0);
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);

    real V = 0;
    #if GEOMETRY == SPHERICAL
      const real x0 = 0.1;
      const real y0 = 0.0;
      const real z0 = 1.0;
    #endif
    for(int k = 0; k < d.np_int[KDIR] ; k++) {
        for(int j = 0; j < d.np_int[JDIR] ; j++) {
            for(int i = 0; i < d.np_int[IDIR] ; i++) {
              #if GEOMETRY == CARTESIAN
                real x = d.x[IDIR](i);
                real y = d.x[JDIR](j);
                real z = d.x[KDIR](k);
                real r=sqrt(x*x+y*y+z*z);

              #elif GEOMETRY == SPHERICAL
                real r = d.x[IDIR](i);
                real th = d.x[JDIR](j);
                real phi = d.x[KDIR](k);

                real x = r*sin(th)*cos(phi) - x0;
                real y = r*sin(th)*sin(phi) - y0;
                real z = r*cos(th) - z0;

                r=sqrt(x*x+y*y+z*z);
              #endif
              if(r<Rstart) {
                V += d.dV(k,j,i);
              }
    }}}

#ifdef WITH_MPI
    MPI_Allreduce(MPI_IN_PLACE, &V, 1, realMPI, MPI_SUM, MPI_COMM_WORLD);
#endif

    real gamma = data.hydro->eos->GetGamma();
    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
              #if GEOMETRY == CARTESIAN
                real x = d.x[IDIR](i);
                real y = d.x[JDIR](j);
                real z = d.x[KDIR](k);
                real r=sqrt(x*x+y*y+z*z);

              #elif GEOMETRY == SPHERICAL
                real r = d.x[IDIR](i);
                real th = d.x[JDIR](j);
                real phi = d.x[KDIR](k);

                real x = r*sin(th)*cos(phi) - x0;
                real y = r*sin(th)*sin(phi) - y0;
                real z = r*cos(th) - z0;

                r=sqrt(x*x+y*y+z*z);
              #endif

              // Sedov Blast Wave Following Stone+2018, 3.4.4
              d.Vc(RHO,k,j,i) = 1.0;
              d.Vc(VX1,k,j,i) = 0.0;
              d.Vc(VX2,k,j,i) = 0.0;
              d.Vc(VX3,k,j,i) = 0.0;

              d.Vc(PRS,k,j,i) = 0.01;

              if(r<Rstart) {
                d.Vc(PRS,k,j,i) = (gamma-1)/V;
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
