// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef ARRAYS_HPP_
#define ARRAYS_HPP_

#include "idefix.hpp"
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
/*
template <typename T> using IdefixHostArray1D = Kokkos::View<T*, Layout, Host>;
template <typename T> using IdefixHostArray2D = Kokkos::View<T**, Layout, Host>;
template <typename T> using IdefixHostArray3D = Kokkos::View<T***, Layout, Host>;
template <typename T> using IdefixHostArray4D = Kokkos::View<T****, Layout, Host>;
*/


#endif // ARRAYS_HPP_
