#ifndef IDEFIX_ARRAYS_HPP_
#define IDEFIX_ARRAYS_HPP_

#include "real_types.h"
#include "idefix.hpp"

template <typename T> using IdefixArray1D = Kokkos::View<T*, Layout, Device>;
template <typename T> using IdefixArray2D = Kokkos::View<T**, Layout, Device>;
template <typename T> using IdefixArray3D = Kokkos::View<T***, Layout, Device>;
template <typename T> using IdefixArray4D = Kokkos::View<T****, Layout, Device>;

template <typename T> using IdefixDualArray1D = Kokkos::DualView<T*, Layout, Device>;
template <typename T> using IdefixDualArray2D = Kokkos::DualView<T**, Layout, Device>;
template <typename T> using IdefixDualArray3D = Kokkos::DualView<T***, Layout, Device>;
template <typename T> using IdefixDualArray4D = Kokkos::DualView<T****, Layout, Device>;



#endif // IDEFIX_ARRAYS


