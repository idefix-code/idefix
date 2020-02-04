#ifndef KOKKOS_TYPES_H_
#define KOKKOS_TYPES_H_

#include <Kokkos_Core.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_View.hpp>

#include <Kokkos_Macros.hpp> // for KOKKOS_ENABLE_XXX

#include <impl/Kokkos_Error.hpp>

#include "shared/real_type.h"

using Device = Kokkos::DefaultExecutionSpace;

enum KokkosLayout {
  KOKKOS_LAYOUT_LEFT,
  KOKKOS_LAYOUT_RIGHT
};

typedef Kokkos::View<real_t****, Device>  DataArray;
typedef DataArray::HostMirror           DataArrayHost;
typedef Kokkos::MDRangePolicy< Kokkos::Rank<4> > mdrange_policy;
