set(idefix_source_file ${idefix_source_file}
  src/hydro/boundary/boundaryloop.hpp
  src/hydro/boundary/hydroboundary.cpp
  src/hydro/boundary/hydroboundary.hpp)

if(Idefix_MPI)
  message(${idefix_source_file})
  set(idefix_source_file ${idefix_source_file}
    src/hydro/boundary/mpi.cpp
    src/hydro/boundary/mpi.hpp
    )
endif()

set(idefix_source_file ${idefix_source_file} PARENT_SCOPE)
