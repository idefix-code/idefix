######################################################################################
# Idefix MHD astrophysical code
#
# Source file src/fluid/RiemannSolver/HDsolvers/CMakeLists.txt
#
# Last modified : 09/2021
#
# Copyright(C) by :
# - Soufiane Baghdadi <soufiane.baghdadi@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2021)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

target_sources(idefix-core
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/hllcHD.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/hllHD.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/roeHD.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/tvdlfHD.hpp
  )
