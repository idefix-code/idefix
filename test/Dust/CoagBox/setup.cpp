#include "idefix.hpp"
#include "setup.hpp"
#include "units.hpp"
#include "coagulation.hpp"

void CheckConservation(DataBlock &data) {
  static bool firstCall = true;
  static real initialDustMass = 0.0;
  int nSpecies = data.dust.size();
  auto dV = data.dV;

#ifdef SINGLE_PRECISION
  const real threshold = 1e-4;
#else
  const real threshold = 1e-13;
#endif

  // Get dust density arrays
  IdefixArray4D<real> dustRho[21];
  for(int n = 0; n < nSpecies; n++) {
    dustRho[n] = data.dust[n]->Vc;
  }

  real dustMass = 0.0;
  idefix_reduce("Total dust mass",
      data.beg[KDIR], data.end[KDIR],
      data.beg[JDIR], data.end[JDIR],
      data.beg[IDIR], data.end[IDIR],
      KOKKOS_LAMBDA (int k, int j, int i, real &m) {

        real rhoDust = 0.0;
        for(int n = 0; n < nSpecies; n++) {
          rhoDust += dustRho[n](RHO,k,j,i);
        }
        m += dV(k,j,i) * rhoDust;
      },
      Kokkos::Sum<real>(dustMass));

#ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &dustMass, 1,
                realMPI, MPI_SUM, MPI_COMM_WORLD);
#endif

  if(firstCall) {
    initialDustMass = dustMass;
    firstCall = false;
    idfx::cout << "Initial total dust mass = " << initialDustMass << std::endl;
  }
  else {
    real err = std::fabs(dustMass - initialDustMass);
    idfx::cout << "Total dust mass = " << dustMass << " | error = " << err << std::endl;
    // Write it in a file 
    // std::ofstream file("conservation.dat", std::ios::app);
    // file << std::setprecision(17) << err << "\n";
    // file.close();   

    if(err > threshold) {
      std::stringstream str;
      str << "Total dust mass is not conserved" << std::endl;
      idfx::cout << std::setprecision(15) << "Initial = " << initialDustMass << " New = " << dustMass << " Error = " << err << std::endl;
      IDEFIX_ERROR(str);
    }
  }

  idfx::cout << "Analysis: done." << std::endl;
}

Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  const bool coag = input.GetOrSet<bool>("Coala","coag",0,0);

  if(coag) {
    const int  kernel  = input.GetOrSet<int>("Coala","kernel",0,0);
    const int  kpol    = input.GetOrSet<int>("Coala","kpol",0,0);
    const real eps     = input.GetOrSet<real>("Coala","eps",0,0)/idfx::units.GetDensity()*idfx::units.GetLength()*idfx::units.GetLength()*idfx::units.GetLength();
    const real coagCFL = input.GetOrSet<real>("Coala","coagCFL",0,0);

    if(kpol!=0) {
      IDEFIX_ERROR("Only kpol = 0 is available when using Coala in Idefix");
    }
    idfx::cout << "Have dust coagulation" << std::endl;
    int nSpecies = data.dust.size();
    idfx::cout << "Evolving nSpecies = " << nSpecies << " dust species" << std::endl;

    CoagulationGlob = new Coagulation(data, kernel, eps, coagCFL);
    data.EnrollUserStepFirst(&CoagSourceTerm);  // Enroll the coagulation source term
  }
}

Setup::~Setup() {
  delete CoagulationGlob;
}