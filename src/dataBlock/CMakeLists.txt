######################################################################################
# Idefix MHD astrophysical code
#
# Source file src/dataBlock/CMakeLists.txt
#
# Last modified : 10/2023
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2021 - 2023)
# - Soufiane Baghdadi <soufiane.baghdadi@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2021)
# - Clément Robert <clement.robert@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2023)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

add_subdirectory(planetarySystem)

target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/coarsen.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/dataBlock.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/dataBlock.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/dataBlockHost.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/dataBlockHost.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/dumpToFile.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/evolveStage.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/fargo.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/fargo.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/makeGeometry.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/stateContainer.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/stateContainer.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/validation.cpp

  )
