#include "idefix.hpp"
#include "setup.hpp"
#include "lookupTable.hpp"
#include "fluid.hpp"
#include "fluid_defs.hpp"
#include "eos.hpp"
#include "units.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "dataBlock.hpp"

void MakeAnalysis(DataBlock &);

namespace {
struct CoolBoxUnits {
  real kB;
  real m_p;
  real rho_unit;
  real vel_unit;
  real len_unit;
};

inline CoolBoxUnits LoadCoolBoxUnits() {
  return {idfx::units.k_B,
          idfx::units.m_p,
          idfx::units.GetDensity(),
          idfx::units.GetVelocity(),
          idfx::units.GetLength()};
}
}  // namespace

// Default constructor
real Tini;
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
  const CoolBoxUnits units = LoadCoolBoxUnits();
  real kB = units.kB;
    real mu = 0.609;
  real m_p = units.m_p;

  real rho_unit = units.rho_unit;
  real vel_unit = units.vel_unit;

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
  // Mirror data on host for analysis reductions.
  DataBlockHost d(data);
  d.SyncFromDevice();

  const CoolBoxUnits units = LoadCoolBoxUnits();

  real kB = units.kB;
  real mu = 0.609;
  real m_p = units.m_p;
  real XH = 0.71;

  real rho_unit = units.rho_unit;
  real vel_unit = units.vel_unit;
  real len_unit = units.len_unit;

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

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = d.beg[IDIR];
  iend = d.end[IDIR];
  jbeg = d.beg[JDIR];
  jend = d.end[JDIR];
  kbeg = d.beg[KDIR];
  kend = d.end[KDIR];

  real numrt = ZERO_F;
  real denom = ZERO_F;
  real vol_tot = ZERO_F;

  for(int k = kbeg; k < kend; k++) {
    for(int j = jbeg; j < jend; j++) {
      for(int i = ibeg; i < iend; i++) {
        real cellVol = d.dV(k,j,i);
        real rho = d.Vc(RHO,k,j,i);
        real prs = d.Vc(PRS,k,j,i);
        real rho_temperature = (mu*m_p/kB) * prs * rho_unit * pow(vel_unit,2);
        numrt += rho_temperature * cellVol;
        denom += rho * rho_unit * cellVol;
        vol_tot += cellVol;
      }
    }
  }

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
