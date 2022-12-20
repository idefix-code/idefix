add_subdirectory(boundary)
add_subdirectory(electromotiveforce)

if(Idefix_MHD)
  add_subdirectory(MHDsolvers)
else()
  add_subdirectory(HDsolvers)
endif()


target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/addNonIdealMHDFlux.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/addSourceTerms.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/axis.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/axis.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/calcCurrent.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/calcParabolicFlux.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/calcRiemannFlux.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/calcRightHandSide.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/checkDivB.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/coarsenFlow.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/convertConsToPrim.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/convertConsToPrimHD.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/convertConsToPrimMHD.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/slopeLimiter.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/shockFlattening.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/shockFlattening.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/evolveStage.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/fluxHD.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/fluxMHD.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/fluid_defs.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/enroll.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/fluid.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/viscosity.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/viscosity.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/thermalDiffusion.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/thermalDiffusion.cpp
  )
