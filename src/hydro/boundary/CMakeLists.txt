set(idefix_source_file ${idefix_source_file}
  ${CMAKE_CURRENT_LIST_DIR}/boundaryloop.hpp
  ${CMAKE_CURRENT_LIST_DIR}/hydroboundary.cpp
  ${CMAKE_CURRENT_LIST_DIR}/hydroboundary.hpp)

if(Idefix_MPI)
  set(idefix_source_file ${idefix_source_file}
    ${CMAKE_CURRENT_LIST_DIR}/mpi.cpp
    ${CMAKE_CURRENT_LIST_DIR}/mpi.hpp
    )
endif()

set(idefix_source_file ${idefix_source_file} PARENT_SCOPE)
