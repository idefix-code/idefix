######################################################################################
# Idefix MHD astrophysical code
#
# Source file src/dataBlock/planetarySystem/CMakeLists.txt
#
# Last modified : 10/2023
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG / CNRS - 2023)
# - and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/planet.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/planet.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/planetarySystem.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/planetarySystem.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/planetStructs.hpp

  )
