add_subdirectory(boundary)
add_subdirectory(braginskii)
add_subdirectory(constrainedTransport)
add_subdirectory(eos)
add_subdirectory(RiemannSolver)
add_subdirectory(tracer)


target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/addNonIdealMHDFlux.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/addSourceTerms.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/calcCurrent.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/calcParabolicFlux.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/calcRightHandSide.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/checkNan.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/checkDivB.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/coarsenFlow.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/convertConsToPrim.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/drag.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/drag.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/evolveStage.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/fluid_defs.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/enroll.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/fluid.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/viscosity.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/viscosity.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/thermalDiffusion.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/thermalDiffusion.cpp
  )
