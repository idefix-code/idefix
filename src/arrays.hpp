// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/arrays.hpp
//
// Last modified : 06/2022
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2020 - 2022)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef ARRAYS_HPP_
#define ARRAYS_HPP_

#include "idefix.hpp"

// =====================================================================================
// Idefix Array Type Definitions
// =====================================================================================

template <typename T> using IdefixArray1D =
                            Kokkos::View<T*, Layout, Device>;
template <typename T> using IdefixArray2D =
                            Kokkos::View<T**, Layout, Device>;
template <typename T> using IdefixArray3D =
                            Kokkos::View<T***, Layout, Device>;
template <typename T> using IdefixArray4D =
                            Kokkos::View<T****, Layout, Device>;

template <typename T> using IdefixHostArray1D =
                            Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace>;
template <typename T> using IdefixHostArray2D =
                            Kokkos::View<T**, Kokkos::LayoutRight, Kokkos::HostSpace>;
template <typename T> using IdefixHostArray3D =
                            Kokkos::View<T***, Kokkos::LayoutRight, Kokkos::HostSpace>;
template <typename T> using IdefixHostArray4D =
                            Kokkos::View<T****, Kokkos::LayoutRight, Kokkos::HostSpace>;

// Atomic arrays
template <typename T> using IdefixAtomicArray1D =
                            Kokkos::View<T*, Layout, Device,
                                         Kokkos::MemoryTraits<Kokkos::Atomic>>;
template <typename T> using IdefixAtomicArray2D =
                            Kokkos::View<T**, Layout, Device,
                                         Kokkos::MemoryTraits<Kokkos::Atomic>>;
template <typename T> using IdefixAtomicArray3D =
                            Kokkos::View<T***, Layout, Device,
                                         Kokkos::MemoryTraits<Kokkos::Atomic>>;

// Note: generic code working with IdefixArrays (Kokkos::View) should use the
// Kokkos traits/members directly (Kokkos::is_view_v, View::rank,
// View::non_const_value_type, View::memory_space,
// Kokkos::SpaceAccessibility<...>::accessible, ...), see Kokkos_View.hpp and
// Kokkos_Concepts.hpp in src/kokkos.

#endif // ARRAYS_HPP_
