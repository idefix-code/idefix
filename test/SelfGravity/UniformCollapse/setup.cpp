#include <cmath>

#include "idefix.hpp"
#include "setup.hpp"

real RcloudGlob;
real internalAreaGlob;

// Analyse data to produce an ascii output
void Analysis(DataBlock & data) {
  // Get the usefull parameters
  real centralMass = data.gravity->centralMass;
  real time = data.t;

  if(time==0.){
    // Write the data in ascii to our file
    std::ofstream f;
    f.open("timevol.dat",std::ios::trunc);
    f.precision(15);
    f << std::scientific << time << "\t" << centralMass << std::endl;
    f.close();
  } else {
    // Append the data in ascii to our file
    std::ofstream f;
    f.open("timevol.dat",std::ios::app);
    f.precision(15);
    f << std::scientific << time << "\t" << centralMass << std::endl;
    f.close();
  }
}

// Compute user variables which will be written in vtk files
void ComputeUserVars(DataBlock & data, UserDefVariablesContainer &variables) {
  // Make references to the user-defined arrays (variables is a container of IdefixHostArray3D)
  // Note that the labels should match the variable names in the input file
  IdefixHostArray3D<real> phiP = variables["phiP"];
  Kokkos::deep_copy(phiP, data.gravity->phiP);
}

// User-defined boundaries
void UserdefBoundary(Hydro *hydro, int dir, BoundarySide side, real t) {
    auto *data = hydro->data;
    if((dir==IDIR) && (side == left)) {
        IdefixArray4D<real> Vc = hydro->Vc;
        int ighost = data->nghost[IDIR];
        IdefixArray1D<real> r = data->x[IDIR];

        idefix_for("UserDefBoundary",
          0, data->np_tot[KDIR],
          0, data->np_tot[JDIR],
          0, ighost,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost);
                        // We cap radial velocity at zero following Xu & Kunz 2021 I.
                        if(Vc(VX1,k,j,ighost)>=ZERO_F) {
                           Vc(VX1,k,j,i) = ZERO_F;
                        } else {
                           Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost);
                        }
                    });
    }
}

void FluxBoundary(Fluid<DefaultPhysics> *hydro, int dir, BoundarySide side, const real t) {
    if((dir==IDIR) && (side == left)) {
      // Loading needed data
      DataBlock &data = *hydro->data;
      IdefixArray4D<real> Flux = hydro->FluxRiemann;
      real halfDt = data.dt/2.; // RK2, dt is actually half at each flux calculation
      int iref = data.nghost[IDIR];
      real rin = data.xbeg[IDIR];
      real internalArea = internalAreaGlob;

      // Storage variables
      real fluxTotMass;

      // Setting positive mass flux to zero to avoid mass flowing in the
      // domain (Xu & Kunz 2021 I.)
      idefix_for("SetZeroOutwardFlux",
                 data.beg[KDIR], data.end[KDIR],
                 data.beg[JDIR], data.end[JDIR],
                 iref, iref+1,
                 KOKKOS_LAMBDA (int k, int j, int i) {
                   if(Flux(RHO,k,j,i)>=ZERO_F) {
                     Flux(RHO,k,j,i)=ZERO_F;
                   }
      });

      idefix_reduce("ComputeTotalMassFlux",
                    data.beg[KDIR], data.end[KDIR],
                    data.beg[JDIR], data.end[JDIR],
                    iref, iref+1,
                    KOKKOS_LAMBDA (int k, int j, int i, real &localSum) {
                      localSum += Flux(RHO,k,j,i);
                    },
                    Kokkos::Sum<real>(fluxTotMass));

      // We normalize by the considered totalA then extrapolate to the whole
      // sphere at given radius (cause we're only 1D)
      data.gravity->centralMass -= fluxTotMass*halfDt*4.*M_PI*rin*rin/internalArea;
    }
}

// Default constructor


// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  data.hydro->EnrollUserDefBoundary(&UserdefBoundary);
  data.hydro->EnrollFluxBoundary(&FluxBoundary);
  output.EnrollUserDefVariables(&ComputeUserVars);
  output.EnrollAnalysis(&Analysis);
  RcloudGlob = 60.;
  internalAreaGlob = grid.xbeg[IDIR]*grid.xbeg[IDIR];
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);

    // Grid and block parameters
    IdefixHostArray4D<real> Vc = d.Vc;
    IdefixHostArray1D<real> r = d.x[IDIR];
    real Rcloud = RcloudGlob;

    // Filling initial profiles
    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                Vc(RHO,k,j,i) = (r(i)<=Rcloud)? 0.001 : 0.00001 ;
                Vc(VX1,k,j,i) = ZERO_F;
            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}
