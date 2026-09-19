######################################################################################
# Idefix MHD astrophysical code
#
# Source file src/gravity/CMakeLists.txt
#
# Last modified : 10/2023
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2021 - 2023)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/gravity.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/gravity.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/laplacian.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/laplacian.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/selfGravity.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/selfGravity.cpp
  )
