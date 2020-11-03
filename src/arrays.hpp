// ********************************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ********************************************************************************************************

#ifndef IDEFIX_ARRAYS_HPP_
#define IDEFIX_ARRAYS_HPP_

#include "idefix.hpp"

template <typename T> using IdefixArray1D = Kokkos::View<T*, Layout, Device>;
template <typename T> using IdefixArray2D = Kokkos::View<T**, Layout, Device>;
template <typename T> using IdefixArray3D = Kokkos::View<T***, Layout, Device>;
template <typename T> using IdefixArray4D = Kokkos::View<T****, Layout, Device>;

template <typename T> using IdefixHostArray1D = Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::HostSpace>;
template <typename T> using IdefixHostArray2D = Kokkos::View<T**, Kokkos::LayoutRight, Kokkos::HostSpace>;
template <typename T> using IdefixHostArray3D = Kokkos::View<T***, Kokkos::LayoutRight, Kokkos::HostSpace>;
template <typename T> using IdefixHostArray4D = Kokkos::View<T****, Kokkos::LayoutRight, Kokkos::HostSpace>;


/*
template <typename T> using IdefixHostArray1D = Kokkos::View<T*, Layout, Host>;
template <typename T> using IdefixHostArray2D = Kokkos::View<T**, Layout, Host>;
template <typename T> using IdefixHostArray3D = Kokkos::View<T***, Layout, Host>;
template <typename T> using IdefixHostArray4D = Kokkos::View<T****, Layout, Host>;
*/


#endif // IDEFIX_ARRAYS


