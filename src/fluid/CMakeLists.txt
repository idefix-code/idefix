######################################################################################
# Idefix MHD astrophysical code
#
# Source file src/fluid/CMakeLists.txt
#
# Last modified : 10/2023
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2023)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

add_subdirectory(boundary)
add_subdirectory(braginskii)
add_subdirectory(constrainedTransport)
add_subdirectory(eos)
add_subdirectory(RiemannSolver)
add_subdirectory(tracer)


target_sources(idefix-core
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/addNonIdealMHDFlux.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/addSourceTerms.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/calcCurrent.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/calcParabolicFlux.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/calcRightHandSide.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/checkNan.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/checkDivB.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/coarsenFlow.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/convertConsToPrim.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/drag.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/drag.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/evolveStage.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/fluid_defs.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/enroll.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/fluid.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/viscosity.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/viscosity.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/thermalDiffusion.hpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/thermalDiffusion.cpp
  )
