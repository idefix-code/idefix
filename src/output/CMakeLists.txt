######################################################################################
# Idefix MHD astrophysical code
#
# Source file src/output/CMakeLists.txt
#
# Last modified : 10/2023
#
# Copyright(C) by :
# - Soufiane Baghdadi <soufiane.baghdadi@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2021)
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2023)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

target_sources(idefix-core
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/slice.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/slice.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/dump.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/dump.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/output.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/output.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/scalarField.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/vtk.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/vtk.hpp
  )
