######################################################################################
# Idefix MHD astrophysical code
#
# Source file src/fluid/RiemannSolver/MHDsolvers/CMakeLists.txt
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
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/hlldMHD.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/hllMHD.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/roeMHD.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/tvdlfMHD.hpp
  )
