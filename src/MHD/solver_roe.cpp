#include "../idefix.hpp"
#include "solvers.hpp"

#define ROE_AVERAGE 0

// Compute Riemann fluxes from states using ROE solver
void Roe(DataBlock & data, int dir, real gamma, real C2Iso) {

    Kokkos::Profiling::pushRegion("ROE_Solver");

    Kokkos::Profiling::popRegion();

}
