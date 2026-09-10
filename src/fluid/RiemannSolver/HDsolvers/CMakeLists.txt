######################################################################################
# Idefix MHD astrophysical code
#
# Source file src/fluid/RiemannSolver/HDsolvers/CMakeLists.txt
#
# Last modified : 09/2021
#
# Copyright(C) by :
# - Soufiane Baghdadi <soufiane.baghdadi@univ-grenoble-alpes.fr> (IPAG / CNRS - 2021)
# - and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/hllcHD.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/hllHD.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/roeHD.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/tvdlfHD.hpp
  )
