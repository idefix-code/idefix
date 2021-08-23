set(idefix_source_file ${idefix_source_file}
  hydro/boundary/boundaryloop.hpp
  hydro/boundary/hydroboundary.cpp
  hydro/boundary/hydroboundary.hpp
  PARENT_SCOPE)

if(MPI)
  set(idefix_source_file ${idefix_source_file}
    hydro/boundary/mpi.cpp
    hydro/boundary/mpi.hpp
    PARENT_SCOPE)
endif()
