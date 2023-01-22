// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************


#ifndef UTILS_VECTOR_HPP_
#define UTILS_VECTOR_HPP_

#include "idefix.hpp"
// a vector structure of size N for reduction in sums
// This struct(=class) defines the += operator on its elements
// Which is used by the reduction function. Other operators
// Can be implemented as well.
template<class T, int N>
struct Vector {
    T v[N];
    KOKKOS_FUNCTION Vector() {
      for (int i = 0; i < N; ++i) {
        v[i] = 0;
      }
    }

    KOKKOS_FUNCTION void operator+=(Vector const volatile& vec) volatile {
        for (int i = 0; i < N; ++i) {
            v[i] = v[i] + vec.v[i];
        }
    }
};

// Shortcut for what follows: a vector of 2 reals
typedef Vector<real,2> MyVector;

// Define the reduction operator in Kokkos space
namespace Kokkos {
template<>
struct reduction_identity< MyVector > {
    KOKKOS_FORCEINLINE_FUNCTION static MyVector sum() {
       return MyVector();
    }
};
}

#endif// UTILS_VECTOR_HPP_
