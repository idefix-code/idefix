// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// These shortcuts are largely inspired from K-Athena
// https://gitlab.com/pgrete/kathena
// by P. Grete.

#ifndef REDUCE_HPP_
#define REDUCE_HPP_

#include <string>
#include "idefix.hpp"
#include "global.hpp"


// 1D default loop pattern
template <typename Function, typename Reducer>
inline void idefix_reduce(const std::string & NAME,
                const int & IB, const int & IE,
                Function function,
                Reducer redFunction) {
    #ifdef DEBUG
    idfx::pushRegion("idefix_reduce("+NAME+")");
    #endif
    Kokkos::parallel_reduce(NAME,
      Kokkos::RangePolicy<>(IB,IE), function, redFunction);
    #ifdef DEBUG
    Kokkos::fence();
    idfx::popRegion();
    #endif
}

// 2D default loop pattern
template <typename Function, typename Reducer>
inline void idefix_reduce(const std::string & NAME,
                const int & JB, const int & JE,
                const int & IB, const int & IE,
                Function function,
                Reducer redFunction) {
    #ifdef DEBUG
    idfx::pushRegion("idefix_reduce("+NAME+")");
    #endif

    // We only implement MDRange reductions here since the other implementations are too
    // complicated to be implemented for any reduction operator on any class
    Kokkos::parallel_reduce(NAME,
      Kokkos::MDRangePolicy<Kokkos::Rank<2, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>
        ({JB,IB},{JE,IE}), function, redFunction);

    #ifdef DEBUG
    Kokkos::fence();
    idfx::popRegion();
    #endif
}

// 3D default loop pattern
template <typename Function, typename Reducer>
inline void idefix_reduce(const std::string & NAME,
                const int & KB, const int & KE,
                const int & JB, const int & JE,
                const int & IB, const int & IE,
                Function function,
                Reducer redFunction) {
    // We only implement MDRange reductions here since the other implementations are too
    // complicated to be implemented for any reduction operator on any class
    #ifdef DEBUG
    idfx::pushRegion("idefix_reduce("+NAME+")");
    #endif
    Kokkos::parallel_reduce(NAME,
      Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>
        ({KB,JB,IB},{KE,JE,IE}), function, redFunction);

    #ifdef DEBUG
    Kokkos::fence();
    idfx::popRegion();
    #endif
}

// 4D default loop pattern
template <typename Function, typename Reducer>
inline void idefix_reduce(const std::string & NAME,
                const int & NB, const int & NE,
                const int & KB, const int & KE,
                const int & JB, const int & JE,
                const int & IB, const int & IE,
                Function function,
                Reducer redFunction) {
    // We only implement MDRange reductions here since the other implementations are too
    // complicated to be implemented for any reduction operator on any class
    #ifdef DEBUG
    idfx::pushRegion("idefix_reduce("+NAME+")");
    #endif
    Kokkos::parallel_reduce(NAME,
      Kokkos::MDRangePolicy<Kokkos::Rank<4, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>
        ({NB,KB,JB,IB},{NE,KE,JE,IE}), function, redFunction);

    #ifdef DEBUG
    Kokkos::fence();
    idfx::popRegion();
    #endif
}

#endif // REDUCE_HPP_
