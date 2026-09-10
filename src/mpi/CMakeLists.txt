######################################################################################
# Idefix MHD astrophysical code
#
# Source file src/mpi/CMakeLists.txt
#
# Last modified : 12/2025
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG / CNRS - 2025)
# - and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/buffer.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/exchanger.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/exchanger.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/mpi.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/mpi.cpp
)
