######################################################################################
# Idefix MHD astrophysical code
#
# Source file src/utils/CMakeLists.txt
#
# Last modified : 01/2025
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG / CNRS - 2022 - 2025)
# - and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

add_subdirectory(iterativesolver)

target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/bigEndian.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/dumpImage.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/dumpImage.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/lookupTable.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/column.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/column.hpp
  )
