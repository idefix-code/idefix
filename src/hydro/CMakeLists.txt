add_subdirectory(boundary)
add_subdirectory(electromotiveforce)
if(Idefix_MHD)
  message("having MHD")
  add_subdirectory(MHDsolvers)
else()
  add_subdirectory(HDsolvers)
endif()


set(idefix_source_file ${idefix_source_file}
  hydro/addNonIdealMHDFlux.hpp
  hydro/addSourceTerms.cpp
  hydro/axis.cpp
  hydro/axis.hpp
  hydro/calcCurrent.cpp
  hydro/calcParabolicFlux.hpp
  hydro/calcRiemannFlux.hpp
  hydro/calcRightHandSide.hpp
  hydro/checkDivB.cpp
  hydro/convertConsToPrim.cpp
  hydro/convertConsToPrimHD.hpp
  hydro/convertConsToPrimMHD.hpp
  hydro/extrapolatePrimVar.hpp
  hydro/fargo.cpp
  hydro/fargo.hpp
  hydro/fluxHD.hpp
  hydro/fluxMHD.hpp
  hydro/hydro_defs.hpp
  hydro/hydro.cpp
  hydro/hydro.hpp
  hydro/viscosity.hpp
  hydro/viscosity.cpp
  PARENT_SCOPE)
