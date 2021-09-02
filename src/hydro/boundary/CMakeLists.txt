target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/boundaryloop.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/hydroboundary.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/hydroboundary.hpp
  )

if(Idefix_MPI)
  target_sources(idefix
    PUBLIC ${CMAKE_CURRENT_LIST_DIR}/mpi.cpp
    PUBLIC ${CMAKE_CURRENT_LIST_DIR}/mpi.hpp
    )
endif()
