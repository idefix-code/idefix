// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/utils/random.hpp
//
// Last modified : 11/2026
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2020 - 2026)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef UTILS_RANDOM_HPP_
#define UTILS_RANDOM_HPP_

#include <cstdint>
#include <limits>

#include "idefix.hpp"
#include "arrays.hpp"

namespace idfx {

/// Self-managed xorshift64 random number generator (see agent/adr/0002).
///
/// One 64-bit state per entry is stored in a flat device array and advanced
/// with the exact same step as Kokkos::Random_XorShift64::urand64(), and the
/// normal draws use the same Box–Muller transform as
/// Kokkos::Random_XorShift64::normal(). The distribution is therefore
/// unchanged, while the sequence becomes deterministic on every backend and
/// thread count and identical on all MPI ranks (unlike
/// Kokkos::Random_XorShift64_Pool, whose no-argument get_state() selects a
/// state through a race on GPU).
///
/// The state is a flat 1D array so that the class stays agnostic of how the
/// caller maps logical indices to streams: a caller with 2D (l, dir) indices
/// simply passes the effective index l*n1+dir (the LayoutRight offset of
/// (l, dir) in a 2D array).
///
/// Because the state lives in a plain Kokkos::View, the class is safe to
/// capture by value in a KOKKOS_LAMBDA: the copy shares the same underlying
/// state, so draws made inside a kernel advance the generator.
///
/// The external draw API is indexed (Drand(i) / Normal(i) on stream i); the
/// raw single-state-word step and Box–Muller transform are private Impl
/// helpers, so the two cannot be confused.
class RandomXorShift64 {
 public:
  RandomXorShift64() = default;
  /// Allocate the state array with n streams. Seeding is done separately
  /// with InitFromSeed().
  explicit RandomXorShift64(int n) : state_("RandomXorShift64State", n) {}

  // The class only holds a View (a reference-counted handle), so the
  // defaulted copy semantics are exactly what device lambdas need.
  RandomXorShift64(const RandomXorShift64&) = default;
  RandomXorShift64& operator=(const RandomXorShift64&) = default;

  /// Deterministically fill the state array from a seed, mirroring the
  /// warm-up and 16-bit packing used by Kokkos::Random_XorShift64_Pool::init
  /// with a single host-side generator, so the initial states are identical
  /// on every rank (host only, defined in random.cpp).
  void InitFromSeed(int seed);

  /// External API: advance state_[i] with the xorshift64 step and return a
  /// double in [0,1), exactly like Kokkos::Random_XorShift64::drand().
  /// Const like Kokkos::View::operator(): the View handle is const but the
  /// underlying state word is advanced in place (this is what makes draws
  /// work inside a by-value KOKKOS_LAMBDA capture).
  KOKKOS_INLINE_FUNCTION
  double Drand(int i) const {
    return DrandImpl(state_(i));
  }

  /// External API: standard normal draw from state_[i] through the
  /// Box–Muller transform, exactly like Kokkos::Random_XorShift64::normal().
  KOKKOS_INLINE_FUNCTION
  double Normal(int i) const {
    return NormalImpl(state_(i));
  }

  /// Raw access to the state array, e.g. to register it in a restart dump.
  IdefixArray1D<uint64_t>& state() { return state_; }
  const IdefixArray1D<uint64_t>& state() const { return state_; }

 private:
  // Internal helpers, operating on a raw 64-bit state word (advanced in
  // place). They are private and carry an Impl suffix so that they cannot be
  // confused with the external (indexed) draw API above.
  /// xorshift64 step, exactly like
  /// Kokkos::Random_XorShift64::urand64()/MAX_URAND64.
  KOKKOS_INLINE_FUNCTION
  static double DrandImpl(uint64_t &state) {
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
    return static_cast<double>((state * 2685821657736338717ULL) - 1)
           / static_cast<double>(Kokkos::Experimental::finite_max_v<uint64_t>);
  }

  /// Box–Muller standard normal draw, exactly like
  /// Kokkos::Random_XorShift64::normal().
  KOKKOS_INLINE_FUNCTION
  static double NormalImpl(uint64_t &state) {
    constexpr double two_pi = 2.0 * Kokkos::numbers::pi_v<double>;
    const double u = DrandImpl(state);
    const double v = DrandImpl(state);
    const double r = Kokkos::sqrt(-2.0 * Kokkos::log(u));
    const double theta = v * two_pi;
    return r * Kokkos::cos(theta);
  }

  IdefixArray1D<uint64_t> state_;
};

}  // namespace idfx

#endif  // UTILS_RANDOM_HPP_
