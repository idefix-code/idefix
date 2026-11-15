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

// Extract dimension from Kokkos::View (1D, 2D, 3D, 4D)
template<typename> struct view_dimension : std::integral_constant<int, 0> {};
template<typename T, typename... Args>
struct view_dimension<Kokkos::View<T*, Args...>> : std::integral_constant<int, 1> {};
template<typename T, typename... Args>
struct view_dimension<Kokkos::View<T**, Args...>> : std::integral_constant<int, 2> {};
template<typename T, typename... Args>
struct view_dimension<Kokkos::View<T***, Args...>> : std::integral_constant<int, 3> {};
template<typename T, typename... Args>
struct view_dimension<Kokkos::View<T****, Args...>> : std::integral_constant<int, 4> {};

// Helper to remove all pointer levels from a type
namespace detail {
template<typename T>
struct remove_all_pointers {
    using type = T;
};

template<typename T>
struct remove_all_pointers<T*> {
    using type = typename remove_all_pointers<T>::type;
};

template<typename T>
struct remove_all_pointers<T* const> {
    using type = typename remove_all_pointers<T>::type;
};

template<typename T>
struct remove_all_pointers<T* volatile> {
    using type = typename remove_all_pointers<T>::type;
};

template<typename T>
struct remove_all_pointers<T* const volatile> {
    using type = typename remove_all_pointers<T>::type;
};  // namespace detail

} // namespace detail

template<typename T>
using remove_all_pointers_t = typename detail::remove_all_pointers<T>::type;

// Extract value_type from Kokkos::View
// For views like Kokkos::View<double***, ...>, we want to extract double, not double***
template<typename> struct view_value_type {};
template<typename T, typename... Args>
struct view_value_type<Kokkos::View<T, Args...>> {
    using type = remove_all_pointers_t<T>;
};

// Check if type is an IdefixArray
template<typename T>
struct is_idefix_array : is_kokkos_view<T> {};

} // namespace detail

// =====================================================================================
// Memory space traits - use detection idiom
// =====================================================================================

namespace detail {

// Helper to check if a type has a specific memory space
template<typename T, typename = void>
struct has_host_space : std::false_type {};

template<typename T>
struct has_host_space<T, std::void_t<typename T::memory_space>> : std::true_type {};

} // namespace detail

// Primary templates
template<typename T>
struct is_host_space : std::false_type {};

template<typename T>
struct is_device_space : std::false_type {};

// Specializations for Kokkos::View
template<typename T, typename... Args>
struct is_host_space<Kokkos::View<T, Args...>> :
    std::disjunction<
        std::is_same<Kokkos::HostSpace, typename Kokkos::View<T, Args...>::memory_space>,
        std::is_same<Kokkos::HostSpace, typename Kokkos::View<T, Args...>::device_type>
    > {};

template<typename T, typename... Args>
struct is_device_space<Kokkos::View<T, Args...>> :
    std::is_same<Device, typename Kokkos::View<T, Args...>::memory_space> {};

// =====================================================================================
// Convenience type traits (public API)
// =====================================================================================

template<typename T>
inline constexpr bool is_kokkos_view_v = idefix::detail::is_kokkos_view<T>::value;
template<typename T>
inline constexpr bool is_idefix_array_v = idefix::detail::is_idefix_array<T>::value;
template<typename T>
inline constexpr int view_dimension_v = idefix::detail::view_dimension<T>::value;
template<typename T>
inline constexpr bool is_host_space_v = idefix::is_host_space<T>::value;
template<typename T>
inline constexpr bool is_device_space_v = idefix::is_device_space<T>::value;
template<typename T>
using value_type_t = typename idefix::detail::view_value_type<T>::type;

} // namespace idefix

#endif // ARRAYS_HPP_
