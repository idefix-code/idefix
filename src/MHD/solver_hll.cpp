#include "../idefix.hpp"
#include "solvers.hpp"

// Compute Riemann fluxes from states using HLL solver
void Hll(DataBlock & data, int dir, real gamma, real C2Iso) {

    Kokkos::Profiling::pushRegion("HLL_Solver");
    

    Kokkos::Profiling::popRegion();

}
