######################################################################################
# Idefix MHD astrophysical code
#
# Source file src/dataBlock/planetarySystem/CMakeLists.txt
#
# Last modified : 10/2023
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2023)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

target_sources(idefix-core
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/planet.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/planet.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/planetarySystem.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/planetarySystem.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/planetStructs.hpp

  )
