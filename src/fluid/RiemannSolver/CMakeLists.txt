
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
