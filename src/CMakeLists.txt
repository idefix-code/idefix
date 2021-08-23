set (CMAKE_CXX_STANDARD 14)

option(Idefix_MHD "enable MHD" OFF)
option(Idefix_MPI "enable Message Passing Interface parallelisation" OFF)

set(idefix_source_file
  arrays.hpp
  error.cpp
  error.hpp
  global.cpp
  global.hpp
  grid.cpp
  grid.hpp
  gridHost.hpp
  gridHost.cpp
  idefix.hpp
  input.cpp
  input.hpp
  loop.hpp
  macros.hpp
  main.cpp
  profiler.cpp
  profiler.hpp
  real_types.hpp
  setup.hpp
  timeIntegrator.hpp
  timeIntegrator.cpp
  )

add_subdirectory(kokkos)
add_subdirectory(hydro)
add_subdirectory(dataBlock)
add_subdirectory(output)
add_subdirectory(rkl)
add_subdirectory("${PROJECT_BINARY_DIR}" "${PROJECT_BINARY_DIR}/local")

add_executable(idefix ${idefix_source_file})
if(Idefix_MHD)
  add_compile_definitions(MHD=YES)
else()
  add_compile_definitions(MHD=NO)
endif()

target_include_directories(idefix PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )
target_include_directories(idefix PUBLIC
                           kokkos/core/src
                           dataBlock/
                           hydro/
                           hydro/boundary
                           hydro/electromotiveforce
                           hydro/HDsolvers
                           hydro/MHDsolvers
                           output/
                           rkl/
                           ./)
target_include_directories(idefix PUBLIC
                          data)
target_include_directories(idefix PUBLIC
                           "${PROJECT_BINARY_DIR}/src/kokkos"
                           )
target_link_libraries(idefix Kokkos::kokkos)
