######################################################################################
# Idefix MHD astrophysical code
#
# Source file src/fluid/braginskii/CMakeLists.txt
#
# Last modified : 10/2023
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2023)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/bragThermalDiffusion.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/bragThermalDiffusion.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/bragViscosity.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/bragViscosity.cpp
  )
