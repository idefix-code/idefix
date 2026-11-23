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

// =====================================================================================
// Type Traits for Kokkos Views / IdefixArrays
// =====================================================================================

namespace idefix {
namespace detail {

// Check if T is a Kokkos::View
template<typename> struct is_kokkos_view : std::false_type {};
template<typename... Args>
struct is_kokkos_view<Kokkos::View<Args...>> : std::true_type {};

// Number of dimensions of a Kokkos::View (0 if T is not a View)
template<typename> struct view_dimension : std::integral_constant<int, 0> {};
template<typename... Args>
struct view_dimension<Kokkos::View<Args...>> :
    std::integral_constant<int, Kokkos::View<Args...>::rank> {};

// Check whether a Kokkos::View lives in host memory
template<typename T>
struct is_host_space : std::false_type {};
template<typename... Args>
struct is_host_space<Kokkos::View<Args...>> :
    std::is_same<Kokkos::HostSpace, typename Kokkos::View<Args...>::memory_space> {};

} // namespace detail

// Value type stored by a Kokkos::View (e.g. double for Kokkos::View<double***>)
template<typename> struct view_value_type {};
template<typename... Args>
struct view_value_type<Kokkos::View<Args...>> {
    using type = typename Kokkos::View<Args...>::non_const_value_type;
};

template<typename T>
inline constexpr bool is_kokkos_view_v = detail::is_kokkos_view<T>::value;
template<typename T>
inline constexpr int view_dimension_v = detail::view_dimension<T>::value;
template<typename T>
inline constexpr bool is_host_space_v = detail::is_host_space<T>::value;
template<typename T>
using value_type_t = typename view_value_type<T>::type;

} // namespace idefix

#endif // ARRAYS_HPP_
