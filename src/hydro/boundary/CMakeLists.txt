set(idefix_source_file ${idefix_source_file}
  src/hydro/boundary/boundaryloop.hpp
  src/hydro/boundary/hydroboundary.cpp
  src/hydro/boundary/hydroboundary.hpp
  PARENT_SCOPE)

if(Idefix_MPI)
  set(idefix_source_file ${idefix_source_file}
    src/hydro/boundary/mpi.cpp
    src/hydro/boundary/mpi.hpp
    PARENT_SCOPE)
endif()
