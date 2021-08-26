cmake_minimum_required(VERSION 3.16)
set(CMAKE_BUILD_TYPE Release)
set (CMAKE_CXX_STANDARD 14)
project (idefix)
option(Idefix_MHD "enable MHD" OFF)
option(Idefix_MPI "enable Message Passing Interface parallelisation" OFF)
set(Idefix_DEFS "definitions.hpp" CACHE STRING "Problem definition header file")
set(idefix_source_file
  ${PROJECT_BINARY_DIR}/setup.cpp)


# load git revision tools
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)

add_subdirectory(src/kokkos)
add_subdirectory(src)

message(STATUS "Idefix Configuration")
# If a CMakeLists.txt is in the problem dir (for problem-specific source files)
# then read it
if(EXISTS ${PROJECT_BINARY_DIR}/CMakeLists.txt)
  add_subdirectory(${PROJECT_BINARY_DIR})
  message(STATUS "Including problem-specific source files in '${PROJECT_BINARY_DIR}'")
endif()

if(Kokkos_ENABLE_CUDA)
  set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "Idefix requires lambdas on Cuda" FORCE)
endif()

add_executable(idefix ${idefix_source_file})


if(Idefix_MHD)
  add_compile_definitions("MHD=YES")
else()
  add_compile_definitions("MHD=NO")
endif()
message(STATUS "    MHD: ${Idefix_MHD}")

if(Idefix_MPI)
  add_compile_definitions("WITH_MPI")
  find_package(MPI REQUIRED)
  target_link_libraries(idefix MPI::MPI_CXX)
endif()
message(STATUS "    MPI: ${Idefix_MPI}")

#update gitversion.hpp if possible
git_describe(GIT_SHA1)
write_file(${CMAKE_SOURCE_DIR}/src/gitversion.hpp "#define GITVERSION \"${GIT_SHA1}\"")
message(STATUS "    Version: ${GIT_SHA1}")

if(NOT ${Idefix_DEFS} STREQUAL "definitions.hpp")
  add_compile_definitions("DEFINITIONS_FILE=\"${Idefix_DEFS}\"")
endif()
message(STATUS "    Problem definitions: '${Idefix_DEFS}'")


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
