// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/utils/random.cpp
//
// Last modified : 11/2026
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2020 - 2026)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "random.hpp"

namespace idfx {

void RandomXorShift64::InitFromSeed(int seed) {
  // Deterministic initialisation of one xorshift64 state per entry from the
  // configured seed (see agent/adr/0002). This mirrors the warm-up and the
  // 16-bit packing used by Kokkos::Random_XorShift64_Pool::init, with a
  // single host-side generator, so the initial states are identical on every
  // rank. A zero state would be stuck at zero, so it is replaced as in the
  // Kokkos generator constructor.
  uint64_t genState = (seed == 0) ? static_cast<uint64_t>(1318319)
                                  : static_cast<uint64_t>(seed);
  auto hostRand = [&genState]() {  // Kokkos::Random_XorShift64::rand() == urand()/2
    genState ^= genState >> 12;
    genState ^= genState << 25;
    genState ^= genState >> 27;
    uint64_t tmp = (genState * 2685821657736338717ULL) >> 16;
    return static_cast<int>((tmp & 0xFFFFFFFFULL) / 2);
  };
  for (int i = 0; i < 17; i++) hostRand();  // warm-up, as in the Kokkos pool
  auto stateHost = Kokkos::create_mirror_view(state_);
  const int n = static_cast<int>(state_.extent(0));
  for (int i = 0; i < n; i++) {
    uint64_t state = 0;
    state |= (static_cast<uint64_t>(hostRand()) & 0xFFFF) << 0;
    state |= (static_cast<uint64_t>(hostRand()) & 0xFFFF) << 16;
    state |= (static_cast<uint64_t>(hostRand()) & 0xFFFF) << 32;
    state |= (static_cast<uint64_t>(hostRand()) & 0xFFFF) << 48;
    stateHost(i) = (state == 0) ? static_cast<uint64_t>(1318319) : state;
  }
  Kokkos::deep_copy(state_, stateHost);
}

}  // namespace idfx
