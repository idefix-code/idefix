if(Idefix_MHD)
  add_subdirectory(MHDsolvers)
else()
  add_subdirectory(HDsolvers)
endif()

target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/calcFlux.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/riemannSolver.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/shockFlattening.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/shockFlattening.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/slopeLimiter.hpp
  )
