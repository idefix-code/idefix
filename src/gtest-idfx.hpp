// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/utils/gtest-idfx.hpp
//
// Last modified : 10/2026
//
// Copyright(C) by :
// - Sébastien Valat <sebastien.valat@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef GTEST_IDFX_HPP_
#define GTEST_IDFX_HPP_

// includes
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <string>
#include "idefix.hpp"

std::string idfx_get_data_file_path(const std::string & sourceFile, const std::string & fname);

/**
 * A macro to load a datafile present at data/{fname}, data/ being in the same dir than the sources.
 * @param fname The file name to append to the SRC_DIR/data path.
*/
#define DATA_FILE(fname) (idfx_get_data_file_path(__FILE__, fname))

/**
 * Compare two C arrays of given size.
 * @param arr1 First C pointer.
 * @param arr2 Second C pointer.
 * @param cnt Number of elements to compare.
 * @param msg Debugging message to print on error.
**/
#define ASSERT_C_ARRAY_EQ(arr1, arr2, cnt, msg) \
  do { \
    for (size_t __i__ = 0 ; __i__ < cnt ; __i__++) \
      ASSERT_EQ(arr1[__i__], arr2[__i__]) << "index=" << __i__ << "\n" << msg; \
  } while(0)

/**
 * Apply a google mock patcher matcher to each element of an array
 * @param array Name of the Kokkos::view to validate.
 * @param matcher Name of the google mock matcher to apply.
 * @param ... Optional arguments to pass to the matcher constructor.
 *
 * For usage example, please see ./src/tests/testKokkos.cpp
 *
 * For what can be used as matcher name, see :
 * https://google.github.io/googletest/reference/matchers.html#generic-comparison
 */
#define ASSERT_THAT_1D_ARRAY(array, matcher, ...) \
  do { \
    auto arrayOnHost = Kokkos::create_mirror_view(array); \
    for (size_t i = 0 ; i < arrayOnHost.extent(0) ; i++) { \
      ASSERT_THAT(arrayOnHost(i), matcher(__VA_ARGS__)) << "      at: (" << i << ")"; \
    } \
  } while(0)

/**
 * Apply a google mock patcher matcher to each element of two array applying
 * the given matcher to each element.
 * @param array Name of the Kokkos::view to validate.
 * @param matcher Name of the google mock matcher to apply. It will recieve the array
 * value as first parameter.
 *
 * For usage example, please see ./src/tests/testKokkos.cpp
 *
 * For what can be used as matcher name, see :
 * https://google.github.io/googletest/reference/matchers.html#generic-comparison
 */
#define ASSERT_THAT_1D_ARRAYS(array, matcher, array2) \
  do { \
    auto arrayOnHost = Kokkos::create_mirror_view(array); \
    auto array2OnHost = Kokkos::create_mirror_view(array2); \
    for (size_t i = 0 ; i < arrayOnHost.extent(0) ; i++) { \
      ASSERT_THAT(arrayOnHost(i), matcher(array2OnHost(i))) << "      at: (" << i << ")"; \
    } \
  } while(0)

/**
 * Apply a google mock patcher matcher to each element of two array applying
 * the given matcher to each element.
 * @param array Name of the Kokkos::view to validate.
 * @param matcher Name of the google mock matcher to apply. It will recieve the array
 * value as first parameter.
 * @param ... Optional extra arguments to pass to the matcher constructor.
 *
 * For usage example, please see ./src/tests/testKokkos.cpp
 *
 * For what can be used as matcher name, see :
 * https://google.github.io/googletest/reference/matchers.html#generic-comparison
 */
#define ASSERT_THAT_1D_ARRAYS_ARGS(array, that, array2, ...) \
  do { \
    auto arrayOnHost = Kokkos::create_mirror_view(array); \
    auto array2OnHost = Kokkos::create_mirror_view(array2); \
    for (size_t i = 0 ; i < arrayOnHost.extent(0) ; i++) { \
      ASSERT_THAT(arrayOnHost(i), that(array2OnHost(i), __VA_ARGS__)) \
        << "      at: (" << i << ")"; \
    } \
  } while(0)

/**
 * Apply a google mock patcher matcher to each element of an array
 * @param array Name of the Kokkos::view to validate.
 * @param matcher Name of the google mock matcher to apply.
 * @param ... Optional arguments to pass to the matcher constructor.
 *
 * For usage example, please see ./src/tests/testKokkos.cpp
 *
 * For what can be used as matcher name, see :
 * https://google.github.io/googletest/reference/matchers.html#generic-comparison
 */
#define ASSERT_THAT_2D_ARRAY(array, matcher, ...) \
  do { \
    auto arrayOnHost = Kokkos::create_mirror_view(array); \
    for (size_t i = 0 ; i < arrayOnHost.extent(0) ; i++) { \
      for (size_t j = 0 ; j < arrayOnHost.extent(1) ; j++) { \
        ASSERT_THAT(arrayOnHost(i, j), matcher(__VA_ARGS__)) \
          << "      at: (" << i << " , " << j << ")"; \
      } \
    } \
  } while(0)

/**
 * Apply a google mock patcher matcher to each element of two array applying
 * the given matcher to each element.
 * @param array Name of the Kokkos::view to validate.
 * @param matcher Name of the google mock matcher to apply. It will recieve the array
 * value as first parameter.
 *
 * For usage example, please see ./src/tests/testKokkos.cpp
 *
 * For what can be used as matcher name, see :
 * https://google.github.io/googletest/reference/matchers.html#generic-comparison
 */
#define ASSERT_THAT_2D_ARRAYS(array, matcher, array2) \
  do { \
    auto arrayOnHost = Kokkos::create_mirror_view(array); \
    auto array2OnHost = Kokkos::create_mirror_view(array2); \
    for (size_t i = 0 ; i < arrayOnHost.extent(0) ; i++) { \
      for (size_t j = 0 ; j < arrayOnHost.extent(1) ; j++) { \
        ASSERT_THAT(arrayOnHost(i, j), matcher(array2OnHost(i, j))) \
          << "      at: (" << i << " , " << j << ")"; \
      } \
    } \
  } while(0)

/**
 * Apply a google mock patcher matcher to each element of two array applying
 * the given matcher to each element.
 * @param array Name of the Kokkos::view to validate.
 * @param matcher Name of the google mock matcher to apply. It will recieve the array
 * value as first parameter.
 * @param ... Optional extra arguments to pass to the matcher constructor.
 *
 * For usage example, please see ./src/tests/testKokkos.cpp
 *
 * For what can be used as matcher name, see :
 * https://google.github.io/googletest/reference/matchers.html#generic-comparison
 */
#define ASSERT_THAT_2D_ARRAYS_ARGS(array, matcher, array2, ...) \
  do { \
    auto arrayOnHost = Kokkos::create_mirror_view(array); \
    auto array2OnHost = Kokkos::create_mirror_view(array2); \
    for (size_t i = 0 ; i < arrayOnHost.extent(0) ; i++) { \
      for (size_t j = 0 ; j < arrayOnHost.extent(1) ; j++) { \
        ASSERT_THAT(arrayOnHost(i, j), matcher(array2OnHost(i, j), __VA_ARGS__)) \
          << "      at: (" << i << " , " << j << ")"; \
      } \
    } \
  } while(0)

/**
 * Apply a google mock patcher matcher to each element of an array
 * @param array Name of the Kokkos::view to validate.
 * @param matcher Name of the google mock matcher to apply.
 * @param ... Optional arguments to pass to the matcher constructor.
 *
 * For usage example, please see ./src/tests/testKokkos.cpp
 *
 * For what can be used as matcher name, see :
 * https://google.github.io/googletest/reference/matchers.html#generic-comparison
 */
#define ASSERT_THAT_3D_ARRAY(array, matcher, ...) \
  do { \
    auto arrayOnHost = Kokkos::create_mirror_view(array); \
    for (size_t i = 0 ; i < arrayOnHost.extent(0) ; i++) { \
      for (size_t j = 0 ; j < arrayOnHost.extent(1) ; j++) { \
        for (size_t k = 0 ; k < arrayOnHost.extent(2) ; k++) { \
          ASSERT_THAT(arrayOnHost(i, j, k), matcher(__VA_ARGS__)) \
            << "      at: (" << i << " , " << j << " , " << k << ")"; \
        } \
      } \
    } \
  } while(0)

/**
 * Apply a google mock patcher matcher to each element of two array applying
 * the given matcher to each element.
 * @param array Name of the Kokkos::view to validate.
 * @param matcher Name of the google mock matcher to apply. It will recieve the array
 * value as first parameter.
 *
 * For usage example, please see ./src/tests/testKokkos.cpp
 *
 * For what can be used as matcher name, see :
 * https://google.github.io/googletest/reference/matchers.html#generic-comparison
 */
#define ASSERT_THAT_3D_ARRAYS(array, matcher, array2) \
  do { \
    auto arrayOnHost = Kokkos::create_mirror_view(array); \
    auto array2OnHost = Kokkos::create_mirror_view(array2); \
    for (size_t i = 0 ; i < arrayOnHost.extent(0) ; i++) { \
      for (size_t j = 0 ; j < arrayOnHost.extent(1) ; j++) { \
        for (size_t k = 0 ; k < arrayOnHost.extent(2) ; k++) { \
          ASSERT_THAT(arrayOnHost(i, j, k), matcher(array2OnHost(i, j, k))) \
            << "      at: (" << i << " , " << j << " , " << k << ")"; \
        } \
      } \
    } \
  } while(0)

/**
 * Apply a google mock patcher matcher to each element of two array applying
 * the given matcher to each element.
 * @param array Name of the Kokkos::view to validate.
 * @param matcher Name of the google mock matcher to apply. It will recieve the array
 * value as first parameter.
 * @param ... Optional extra arguments to pass to the matcher constructor.
 *
 * For usage example, please see ./src/tests/testKokkos.cpp
 *
 * For what can be used as matcher name, see :
 * https://google.github.io/googletest/reference/matchers.html#generic-comparison
 */
#define ASSERT_THAT_3D_ARRAYS_ARGS(array, matcher, array2, ...) \
  do { \
    auto arrayOnHost = Kokkos::create_mirror_view(array); \
    auto array2OnHost = Kokkos::create_mirror_view(array2); \
    for (size_t i = 0 ; i < arrayOnHost.extent(0) ; i++) { \
      for (size_t j = 0 ; j < arrayOnHost.extent(1) ; j++) { \
        for (size_t k = 0 ; k < arrayOnHost.extent(2) ; k++) { \
          ASSERT_THAT(arrayOnHost(i, j, k), matcher(array2OnHost(i, j, k), __VA_ARGS__)) \
            << "      at: (" << i << " , " << j << " , " << k << ")"; \
        } \
      } \
    } \
  } while(0)

/**
 * Apply a google mock patcher matcher to each element of an array
 * @param array Name of the Kokkos::view to validate.
 * @param matcher Name of the google mock matcher to apply.
 * @param ... Optional arguments to pass to the matcher constructor.
 *
 * For usage example, please see ./src/tests/testKokkos.cpp
 *
 * For what can be used as matcher name, see :
 * https://google.github.io/googletest/reference/matchers.html#generic-comparison
 */
#define ASSERT_THAT_4D_ARRAY(array, matcher, ...) \
  do { \
    auto arrayOnHost = Kokkos::create_mirror_view(array); \
    for (size_t i = 0 ; i < arrayOnHost.extent(0) ; i++) { \
      for (size_t j = 0 ; j < arrayOnHost.extent(1) ; j++) { \
        for (size_t k = 0 ; k < arrayOnHost.extent(2) ; k++) { \
          for (size_t l = 0 ; l < arrayOnHost.extent(3) ; l++) { \
            ASSERT_THAT(arrayOnHost(i, j, k, l), matcher(__VA_ARGS__)) \
              << "      at: (" << i << " , " << j << " , " << k << " , " << l << ")"; \
          } \
        } \
      } \
    } \
  } while(0)

/**
 * Apply a google mock patcher matcher to each element of two array applying
 * the given matcher to each element.
 * @param array Name of the Kokkos::view to validate.
 * @param matcher Name of the google mock matcher to apply. It will recieve the array
 * value as first parameter.
 *
 * For usage example, please see ./src/tests/testKokkos.cpp
 *
 * For what can be used as matcher name, see :
 * https://google.github.io/googletest/reference/matchers.html#generic-comparison
 */
#define ASSERT_THAT_4D_ARRAYS(array, matcher, array2) \
  do { \
    auto arrayOnHost = Kokkos::create_mirror_view(array); \
    auto array2OnHost = Kokkos::create_mirror_view(array2); \
    for (size_t i = 0 ; i < arrayOnHost.extent(0) ; i++) { \
      for (size_t j = 0 ; j < arrayOnHost.extent(1) ; j++) { \
        for (size_t k = 0 ; k < arrayOnHost.extent(2) ; k++) { \
          for (size_t l = 0 ; l < arrayOnHost.extent(3) ; l++) { \
            ASSERT_THAT(arrayOnHost(i, j, k, l), matcher(array2OnHost(i, j, k, l))) \
              << "      at: (" << i << " , " << j << " , " << k << " , " << l << ")"; \
          } \
        } \
      } \
    } \
  } while(0)

/**
 * Apply a google mock patcher matcher to each element of two array applying
 * the given matcher to each element.
 * @param array Name of the Kokkos::view to validate.
 * @param matcher Name of the google mock matcher to apply. It will recieve the array
 * value as first parameter.
 * @param ... Optional extra arguments to pass to the matcher constructor.
 *
 * For usage example, please see ./src/tests/testKokkos.cpp
 *
 * For what can be used as matcher name, see :
 * https://google.github.io/googletest/reference/matchers.html#generic-comparison
 */
#define ASSERT_THAT_4D_ARRAYS_ARGS(array, matcher, array2, ...) \
  do { \
    auto arrayOnHost = Kokkos::create_mirror_view(array); \
    auto array2OnHost = Kokkos::create_mirror_view(array2); \
    for (size_t i = 0 ; i < arrayOnHost.extent(0) ; i++) { \
      for (size_t j = 0 ; j < arrayOnHost.extent(1) ; j++) { \
        for (size_t k = 0 ; k < arrayOnHost.extent(2) ; k++) { \
          for (size_t l = 0 ; l < arrayOnHost.extent(3) ; l++) { \
            ASSERT_THAT(arrayOnHost(i, j, k, l), matcher(array2OnHost(i, j, k, l), __VA_ARGS__)) \
              << "      at: (" << i << " , " << j << " , " << k << " , " << l << ")"; \
          } \
        } \
      } \
    } \
  } while(0)

#endif //GTEST_IDFX_HPP_
