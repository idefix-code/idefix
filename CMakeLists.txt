cmake_minimum_required(VERSION 3.10)
set(CMAKE_BUILD_TYPE Release)
set (CMAKE_CXX_STANDARD 14)
project (idefix)
option(Idefix_MHD "enable MHD" OFF)
option(Idefix_MPI "enable Message Passing Interface parallelisation" OFF)
set(idefix_source_file
  ${PROJECT_BINARY_DIR}/setup.cpp)


add_subdirectory(src/kokkos)
add_subdirectory(src)


if(Kokkos_ENABLE_CUDA)
  set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "Idefix requires lambdas on Cuda" FORCE)
endif()


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
                           src/kokkos/core/src
                           src/dataBlock/
                           src/hydro/
                           src/hydro/boundary
                           src/hydro/electromotiveforce
                           src/hydro/HDsolvers
                           src/hydro/MHDsolvers
                           src/output/
                           src/rkl/
                           src/)

target_link_libraries(idefix Kokkos::kokkos)
