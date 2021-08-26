add_subdirectory(boundary)
add_subdirectory(electromotiveforce)
if(Idefix_MHD)
  add_subdirectory(MHDsolvers)
else()
  add_subdirectory(HDsolvers)
endif()


set(idefix_source_file ${idefix_source_file}
  src/hydro/addNonIdealMHDFlux.hpp
  src/hydro/addSourceTerms.cpp
  src/hydro/axis.cpp
  src/hydro/axis.hpp
  src/hydro/calcCurrent.cpp
  src/hydro/calcParabolicFlux.hpp
  src/hydro/calcRiemannFlux.hpp
  src/hydro/calcRightHandSide.hpp
  src/hydro/checkDivB.cpp
  src/hydro/convertConsToPrim.cpp
  src/hydro/convertConsToPrimHD.hpp
  src/hydro/convertConsToPrimMHD.hpp
  src/hydro/extrapolatePrimVar.hpp
  src/hydro/fargo.cpp
  src/hydro/fargo.hpp
  src/hydro/fluxHD.hpp
  src/hydro/fluxMHD.hpp
  src/hydro/hydro_defs.hpp
  src/hydro/hydro.cpp
  src/hydro/hydro.hpp
  src/hydro/viscosity.hpp
  src/hydro/viscosity.cpp
  PARENT_SCOPE)
