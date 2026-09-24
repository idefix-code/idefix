######################################################################################
# Idefix MHD astrophysical code
#
# Source file src/fluid/RiemannSolver/CMakeLists.txt
#
# Last modified : 10/2023
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2023)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

add_subdirectory(MHDsolvers)
add_subdirectory(HDsolvers)
add_subdirectory(Dustsolvers)


target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/calcFlux.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/extrapolateToFaces.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/flux.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/riemannSolver.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/shockFlattening.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/slopeLimiter.hpp
  )
