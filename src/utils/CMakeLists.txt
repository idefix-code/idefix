######################################################################################
# Idefix MHD astrophysical code
#
# Source file src/utils/CMakeLists.txt
#
# Last modified : 01/2025
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2022 - 2025)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

add_subdirectory(iterativesolver)

target_sources(idefix-core
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/bigEndian.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/dumpImage.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/dumpImage.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/lookupTable.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/column.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/column.hpp
  )

if (Idefix_GTEST)
  target_sources(unittests
    PRIVATE ${CMAKE_CURRENT_LIST_DIR}/tests/testLookupTable.cpp
  )
endif()
