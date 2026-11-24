// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/utils/tests/testColumn.hpp
//
// Last modified : 10/2026
//
// Copyright(C) by :
// - Sébastien Valat <sebastien.valat@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef UTILS_ARRAYHELPERS_HPP_
#define UTILS_ARRAYHELPERS_HPP_

// headers
#include <type_traits>
#include "idefix.hpp"
#include "arrays.hpp"
#include "gtest-idfx.hpp"
#include "unitTestApp.hpp"
#include "column.hpp"

/**
 * @brief Uniformaly fill the given array with the given value.
 *
 * Uniformaly fill the given array with the given value.
 *
 * @param array The array to fill, being a 1D, 2D, 3D or 3D Kokkos array.
 * @param value The value to assign to each cells.
 */
template <class T, class U>
void ArrayUniformFill(T & array, const U value) {
  if constexpr (array.rank() == 1) {
    idefix_for("ArrayUniformFill.1D", 0, array.extent(0),
      KOKKOS_LAMBDA(int i) {
        array(i) = value;
      }
    );
  } else if constexpr (array.rank() == 2) {
    idefix_for("ArrayUniformFill.2D",
      0, array.extent(0),
      0, array.extent(1),
      KOKKOS_LAMBDA(int j, int i) {
        array(j, i) = value;
      }
    );
  } else if constexpr (array.rank() == 3) {
    idefix_for("ArrayUniformFill.3D",
      0, array.extent(0),
      0, array.extent(1),
      0, array.extent(2),
      KOKKOS_LAMBDA(int k, int j, int i) {
        array(k, j, i) = value;
      }
    );
  } else if constexpr (array.rank() == 4) {
    idefix_for("ArrayUniformFill.3D",
      0, array.extent(0),
      0, array.extent(1),
      0, array.extent(2),
      0, array.extent(3),
      KOKKOS_LAMBDA(int l, int k, int j, int i) {
        array(l, k, j, i) = value;
      }
    );
  } else {
    IDEFIX_ERROR("Unsupported number of dimensions !");
  }
}

#endif //UTILS_ARRAYHELPERS_HPP_
