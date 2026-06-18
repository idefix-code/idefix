#include "idefix.hpp"
#include "setup.hpp"
#include "lookupTable.hpp"
#include "fluid.hpp"
#include "fluid_defs.hpp"
#include "eos.hpp"
#include "input.hpp"
#include "grid.hpp"

void MakeAnalysis(DataBlock &);

// Default constructor
real Tini;
#include "units.hpp"
EquationOfState eos;


// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  Tini    = input.Get<real>("Setup", "Tini", 0);
  eos =  EquationOfState(input, &data, "Hydro");
  output.EnrollAnalysis(&MakeAnalysis);
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    real kB = idfx::units.k_B;
    real mu = 0.609;
    real m_p = idfx::units.m_p;
    real XH = 0.71;

    real rho_unit = idfx::units.GetDensity();
    real vel_unit = idfx::units.GetVelocity();

    // Create a host copy
    DataBlockHost d(data);


    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                d.Vc(RHO,k,j,i)     = ONE_F;
                EXPAND(\
                d.Vc(VX1,k,j,i) =  ONE_F; ,\
                d.Vc(VX2,k,j,i) = ZERO_F;,\
                d.Vc(VX3,k,j,i) = ZERO_F;)
                d.Vc(TRG,k,j,i) = ONE_F;
#if HAVE_ENERGY
                d.Vc(PRS,k,j,i) = ((d.Vc(RHO,k,j,i)*rho_unit)/(mu*m_p) * kB * Tini)/(rho_unit*pow(vel_unit,2));
#endif
            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {
  // Mirror data on Host
  DataBlockHost d(data);

  // Sync it
  d.SyncFromDevice();

  real kB = idfx::units.k_B;
  real mu = 0.609;
  real m_p = idfx::units.m_p;
  real XH = 0.71;

  real rho_unit = idfx::units.GetDensity();
  real vel_unit = idfx::units.GetVelocity();
  real len_unit = idfx::units.GetLength();

  // Ensure time considered is unique
  static real time_old = -ONE_F;
  static real e_now = eos.GetInternalEnergy((rho_unit/(mu*m_p))*kB*Tini, rho_unit); // cgs
  static real e_old = eos.GetInternalEnergy((rho_unit/(mu*m_p))*kB*Tini, rho_unit); // cgs
  real time_now = data.t;
  std::ios_base::openmode mode = std::ios::out | ((time_old==-1.0)?std::ios::trunc : std::ios::app);
  real dt = (time_old>0)?((time_now-time_old)*(len_unit/vel_unit)):ZERO_F;
  if (time_now==time_old) {
    return;
  }
  else {
    time_old = time_now;
    e_old = e_now;
  }


  IdefixArray1D<real> x1 = d.x[IDIR];
  IdefixArray1D<real> x2 = d.x[JDIR];
  IdefixArray1D<real> x3 = d.x[KDIR];
  IdefixArray3D<real> dV = d.dV;

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = d.beg[IDIR];
  iend = d.end[IDIR];
  jbeg = d.beg[JDIR];
  jend = d.end[JDIR];
  kbeg = d.beg[KDIR];
  kend = d.end[KDIR];
/*
  IdefixArray3D<real> rho_temperature("rhotemperatureArray", iend-ibeg+1, jend-jbeg+1, kend-kbeg+1);
  idefix_for("rho_temperatureLoop",
           kbeg,kend,
           jbeg,jend,
           ibeg,iend,
           KOKKOS_LAMBDA (int k, int j, int i) {
              rho_temperature(k,j,i) = (mu*mp/kB)*d.Vc(PRS,k,j,i)*rho_unit*pow(vel_unit,2);
            });
*/
  real numrt = ZERO_F;
  real denom = ZERO_F;
  real vol_tot = ZERO_F;

  idefix_reduce("numerator",
              kbeg,kend,
              jbeg,jend,
              ibeg,iend,
              KOKKOS_LAMBDA (int k, int j, int i, real &numrt) {
                  real rho_temperature = (mu*m_p/kB)*d.Vc(PRS,k,j,i)*rho_unit*pow(vel_unit,2);
                  numrt += (rho_temperature*dV(k,j,i));
              },
              Kokkos::Sum<real> (numrt));
  idefix_reduce("denomenator",
              kbeg,kend,
              jbeg,jend,
              ibeg,iend,
              KOKKOS_LAMBDA (int k, int j, int i, real &denom) {
                  denom += (d.Vc(RHO,k,j,i)*rho_unit*dV(k,j,i));
              },
              Kokkos::Sum<real> (denom));
    idefix_reduce("volume",
              kbeg,kend,
              jbeg,jend,
              ibeg,iend,
              KOKKOS_LAMBDA (int k, int j, int i, real &vol_tot) {
                  vol_tot += dV(k,j,i);
              },
              Kokkos::Sum<real> (vol_tot));

  real temperature_mass_avg = numrt/denom; // K
  real dens_avg = denom/vol_tot; // cgs
  e_now = eos.GetInternalEnergy((dens_avg/(mu*m_p))*kB*temperature_mass_avg, dens_avg); // cgs
  real de = e_now-e_old;

  real de_dt = (dt>0)?(de/dt):ZERO_F;
  real nH = dens_avg*XH/m_p;

  // Write the data in ascii to our file
  std::ofstream file;
  file.open("analysis-time-temperature.dat",mode);
  if (!file.is_open()) {
    IDEFIX_ERROR ("Problem opening analysis file!");
  }
  file.precision(10);
  file << std::scientific << data.t << '\t' << temperature_mass_avg << '\t' << (-de_dt/pow(nH,2)) << std::endl;
  file.close();
}
