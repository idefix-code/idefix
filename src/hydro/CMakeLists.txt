add_subdirectory(boundary)
add_subdirectory(electromotiveforce)
if(Idefix_MHD)
  add_subdirectory(MHDsolvers)
else()
  add_subdirectory(HDsolvers)
endif()


set(idefix_source_file ${idefix_source_file}
  ${CMAKE_CURRENT_LIST_DIR}/addNonIdealMHDFlux.hpp
  ${CMAKE_CURRENT_LIST_DIR}/addSourceTerms.cpp
  ${CMAKE_CURRENT_LIST_DIR}/axis.cpp
  ${CMAKE_CURRENT_LIST_DIR}/axis.hpp
  ${CMAKE_CURRENT_LIST_DIR}/calcCurrent.cpp
  ${CMAKE_CURRENT_LIST_DIR}/calcParabolicFlux.hpp
  ${CMAKE_CURRENT_LIST_DIR}/calcRiemannFlux.hpp
  ${CMAKE_CURRENT_LIST_DIR}/calcRightHandSide.hpp
  ${CMAKE_CURRENT_LIST_DIR}/checkDivB.cpp
  ${CMAKE_CURRENT_LIST_DIR}/convertConsToPrim.cpp
  ${CMAKE_CURRENT_LIST_DIR}/convertConsToPrimHD.hpp
  ${CMAKE_CURRENT_LIST_DIR}/convertConsToPrimMHD.hpp
  ${CMAKE_CURRENT_LIST_DIR}/extrapolatePrimVar.hpp
  ${CMAKE_CURRENT_LIST_DIR}/fargo.cpp
  ${CMAKE_CURRENT_LIST_DIR}/fargo.hpp
  ${CMAKE_CURRENT_LIST_DIR}/fluxHD.hpp
  ${CMAKE_CURRENT_LIST_DIR}/fluxMHD.hpp
  ${CMAKE_CURRENT_LIST_DIR}/hydro_defs.hpp
  ${CMAKE_CURRENT_LIST_DIR}/hydro.cpp
  ${CMAKE_CURRENT_LIST_DIR}/hydro.hpp
  ${CMAKE_CURRENT_LIST_DIR}/viscosity.hpp
  ${CMAKE_CURRENT_LIST_DIR}/viscosity.cpp
  PARENT_SCOPE)
