cmake_minimum_required(VERSION 3.16)
set(CMAKE_BUILD_TYPE Release)
set (CMAKE_CXX_STANDARD 14)

project (idefix VERSION 1.0.0)
option(Idefix_MHD "enable MHD" OFF)
option(Idefix_MPI "enable Message Passing Interface parallelisation" OFF)
option(Idefix_HIGH_ORDER_FARGO "Force Fargo to use a PPM reconstruction scheme" OFF)
set(Idefix_CXX_FLAGS "" CACHE STRING "Additional compiler/linker flag")
set(Idefix_DEFS "definitions.hpp" CACHE STRING "Problem definition header file")


# load git revision tools
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)
include(ReplaceIdefixSource)
include(AddIdefixSource)

#Idefix requires Cuda Lambdas (experimental)
if(Kokkos_ENABLE_CUDA)
  set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "Idefix requires lambdas on Cuda" FORCE)
endif()

# Add kokkos CMAKE files (required early since these set compiler options)
add_subdirectory(src/kokkos build/kokkos)

# Add Idefix CXX Flags
add_compile_options(${Idefix_CXX_FLAGS})

# sub directories should be called after add_executable
# because the binary target "idefix" should be defined
# before adding source files via "target_sources"
#create target binary "idefix"
add_executable(idefix ${PROJECT_BINARY_DIR}/setup.cpp)
add_subdirectory(src build)

# If a CMakeLists.txt is in the problem dir (for problem-specific source files)
# then read it
if(EXISTS ${PROJECT_BINARY_DIR}/CMakeLists.txt)
  message(STATUS "Including problem-specific CMakeLists in '${PROJECT_BINARY_DIR}'")
  add_subdirectory(${PROJECT_BINARY_DIR} build/setup)
endif()



if(Idefix_MHD)
  add_compile_definitions("MHD=YES")
else()
  add_compile_definitions("MHD=NO")
endif()

if(Idefix_MPI)
  add_compile_definitions("WITH_MPI")
  find_package(MPI REQUIRED)
  target_link_libraries(idefix MPI::MPI_CXX)
endif()

if(Idefix_HIGH_ORDER_FARGO)
  add_compile_definitions("HIGH_ORDER_FARGO")
endif()

#update gitversion.hpp if possible
git_describe(GIT_SHA1)
write_file(${CMAKE_SOURCE_DIR}/src/gitversion.hpp "#define GITVERSION \"${GIT_SHA1}\"")

if(NOT ${Idefix_DEFS} STREQUAL "definitions.hpp")
  add_compile_definitions("DEFINITIONS_FILE=\"${Idefix_DEFS}\"")
endif()


target_include_directories(idefix PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )
target_include_directories(idefix PUBLIC
                           src/kokkos/core/src
                           src/dataBlock
                           src/hydro
                           src/hydro/boundary
                           src/hydro/electromotiveforce
                           src/hydro/HDsolvers
                           src/hydro/MHDsolvers
                           src/output
                           src/rkl
                           src
                           )

target_link_libraries(idefix Kokkos::kokkos)

message(STATUS "Idefix final configuration")
message(STATUS "    MHD: ${Idefix_MHD}")
message(STATUS "    MPI: ${Idefix_MPI}")
message(STATUS "    Version: ${GIT_SHA1}")
message(STATUS "    Problem definitions: '${Idefix_DEFS}'")
