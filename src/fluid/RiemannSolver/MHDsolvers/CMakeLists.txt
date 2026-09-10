######################################################################################
# Idefix MHD astrophysical code
#
# Source file src/fluid/RiemannSolver/MHDsolvers/CMakeLists.txt
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
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/hlldMHD.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/hllMHD.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/roeMHD.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/tvdlfMHD.hpp
  )
