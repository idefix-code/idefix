######################################################################################
# Idefix MHD astrophysical code
#
# Source file src/CMakeLists.txt
#
# Last modified : 04/2025
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2021 - 2025)
# - Soufiane Baghdadi <soufiane.baghdadi@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2021)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

add_subdirectory(fluid)
add_subdirectory(dataBlock)
add_subdirectory(output)
add_subdirectory(rkl)
add_subdirectory(gravity)
add_subdirectory(utils)

target_sources(idefix
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/main.cpp
)

target_sources(idefix-core
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/arrays.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/error.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/error.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/global.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/global.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/grid.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/grid.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/gridHost.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/gridHost.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/idefix.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/input.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/input.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/loop.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/macros.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/profiler.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/profiler.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/real_types.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/reduce.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/setup.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/setup.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/units.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/units.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/timeIntegrator.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/timeIntegrator.cpp
)
# NOTE: mpi file are added by the root CMakeLists along with dependencies to the MPI library

if (Idefix_GTEST)
  target_sources(unittests
    PRIVATE ${CMAKE_CURRENT_LIST_DIR}/tests/testKokkos.cpp
    PRIVATE ${CMAKE_CURRENT_LIST_DIR}/tests/testArrays.cpp
    PRIVATE ${CMAKE_CURRENT_LIST_DIR}/gtest-idfx.cpp
  )
endif()

#add_subdirectory("${PROJECT_BINARY_DIR}" "${PROJECT_BINARY_DIR}/local")
