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

target_sources(idefix-core
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/gravity.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/gravity.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/laplacian.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/laplacian.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/selfGravity.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/selfGravity.cpp
  )
