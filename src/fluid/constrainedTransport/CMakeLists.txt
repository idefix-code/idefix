######################################################################################
# Idefix MHD astrophysical code
#
# Source file src/fluid/constrainedTransport/CMakeLists.txt
#
# Last modified : 12/2025
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2023 - 2025)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

target_sources(idefix-core
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/calcCornerEmf.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/calcNonidealEMF.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/calcRiemannEmf.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/constrainedTransport.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/EMFexchange.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/enforceEMFBoundary.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/enforceVectorPotentialBoundary.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/evolveMagField.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/evolveVectorPotential.hpp
  )
