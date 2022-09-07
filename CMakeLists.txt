cmake_minimum_required(VERSION 3.16)
set(CMAKE_BUILD_TYPE Release)
set (CMAKE_CXX_STANDARD 17)

project (idefix VERSION 1.0.0)
option(Idefix_MHD "enable MHD" OFF)
option(Idefix_MPI "enable Message Passing Interface parallelisation" OFF)
option(Idefix_HIGH_ORDER_FARGO "Force Fargo to use a PPM reconstruction scheme" OFF)
option(Idefix_DEBUG "Enable Idefix debug features (makes the code very slow)" OFF)
option(Idefix_WERROR "Treat compiler warnings as errors" OFF)
set(Idefix_CXX_FLAGS "" CACHE STRING "Additional compiler/linker flag")
set(Idefix_DEFS "definitions.hpp" CACHE FILEPATH "Problem definition header file")
set(Idefix_RECONSTRUCTION "Linear" CACHE STRING "Type of cell reconstruction scheme")
if(Idefix_MHD)
  option(Idefix_EVOLVE_VECTOR_POTENTIAL "Evolve the vector potential instead of the field (helps reducing div(B) in long runs)" OFF)
endif()
set_property(CACHE Idefix_RECONSTRUCTION PROPERTY STRINGS Constant Linear LimO3 Parabolic)
set(Idefix_PRECISION "Double" CACHE STRING "Precision of arithmetics")
set_property(CACHE Idefix_PRECISION PROPERTY STRINGS Double Single)


# load git revision tools
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)
include(ReplaceIdefixSource)
include(AddIdefixSource)
include(SetIdefixProperty)

#Idefix requires Cuda Lambdas (experimental)
if(Kokkos_ENABLE_CUDA)
  set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "Idefix requires lambdas on Cuda" FORCE)
endif()

# Add kokkos CMAKE files (required early since these set compiler options)
add_subdirectory(src/kokkos build/kokkos)

# Add Idefix CXX Flags
add_compile_options(${Idefix_CXX_FLAGS})

if(Idefix_WERROR)
  if(Kokkos_ENABLE_CUDA)
    add_compile_options(-Xcudafe --promote_warnings)
  else()
    add_compile_options(-Werror)
  endif()
endif()

# sub directories should be called after add_executable
# because the binary target "idefix" should be defined
# before adding source files via "target_sources"
# create target binary "idefix"
add_executable(idefix)
add_subdirectory(src build)

if(EXISTS ${PROJECT_BINARY_DIR}/setup.cpp)
  target_sources(idefix PUBLIC ${PROJECT_BINARY_DIR}/setup.cpp)
else()
  message(WARNING "No specific setup.cpp found in the problem directory")
  target_sources(idefix PUBLIC ${CMAKE_CURRENT_LIST_DIR}/src/setup.cpp)
endif()

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

if(Idefix_DEBUG)
  add_compile_definitions("DEBUG")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
endif()

if(Idefix_HIGH_ORDER_FARGO)
  add_compile_definitions("HIGH_ORDER_FARGO")
endif()

if(Idefix_EVOLVE_VECTOR_POTENTIAL)
  add_compile_definitions("EVOLVE_VECTOR_POTENTIAL")
endif()
#update gitversion.hpp if possible
git_describe(GIT_SHA1)
write_file(${CMAKE_SOURCE_DIR}/src/gitversion.hpp "#define GITVERSION \"${GIT_SHA1}\"")

if(NOT ${Idefix_DEFS} STREQUAL "definitions.hpp")
  add_compile_definitions("DEFINITIONS_FILE=\"${Idefix_DEFS}\"")
endif()

# Order of the scheme
if(${Idefix_RECONSTRUCTION} STREQUAL "Constant")
  add_compile_definitions("ORDER=1")
elseif(${Idefix_RECONSTRUCTION} STREQUAL "Linear")
  add_compile_definitions("ORDER=2")
elseif(${Idefix_RECONSTRUCTION} STREQUAL "LimO3")
  add_compile_definitions("ORDER=3")
elseif(${Idefix_RECONSTRUCTION} STREQUAL "Parabolic")
  add_compile_definitions("ORDER=4")
else()
  message(ERROR "Reconstruction type '${Idefix_RECONSTRUCTION}' is invalid")
endif()

if(${Idefix_PRECISION} STREQUAL "Single")
  add_compile_definitions("SINGLE_PRECISION")
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
                           src/gravity
                           src/utils
                           src
                           )

target_link_libraries(idefix Kokkos::kokkos)

message(STATUS "Idefix final configuration")
if(Idefix_EVOLVE_VECTOR_POTENTIAL)
  message(STATUS "    MHD: ${Idefix_MHD} (Vector potential)")
else()
  message(STATUS "    MHD: ${Idefix_MHD}")
endif()
message(STATUS "    MPI: ${Idefix_MPI}")
message(STATUS "    Reconstruction: ${Idefix_RECONSTRUCTION}")
message(STATUS "    Precision: ${Idefix_PRECISION}")
message(STATUS "    Version: ${GIT_SHA1}")
message(STATUS "    Problem definitions: '${Idefix_DEFS}'")
