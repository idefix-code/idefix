cmake_minimum_required(VERSION 3.16)
set(CMAKE_BUILD_TYPE Release)
set (CMAKE_CXX_STANDARD 17)

set(Idefix_VERSION_MAJOR 2)
set(Idefix_VERSION_MINOR 1)
set(Idefix_VERSION_PATCH 00)

project (idefix VERSION 2.1.00)
option(Idefix_MHD "enable MHD" OFF)
option(Idefix_MPI "enable Message Passing Interface parallelisation" OFF)
option(Idefix_HIGH_ORDER_FARGO "Force Fargo to use a PPM reconstruction scheme" OFF)
option(Idefix_DEBUG "Enable Idefix debug features (makes the code very slow)" OFF)
option(Idefix_RUNTIME_CHECKS "Enable runtime sanity checks" OFF)
option(Idefix_WERROR "Treat compiler warnings as errors" OFF)
set(Idefix_CXX_FLAGS "" CACHE STRING "Additional compiler/linker flag")
set(Idefix_DEFS "definitions.hpp" CACHE FILEPATH "Problem definition header file")
option(Idefix_CUSTOM_EOS "Use custom equation of state" OFF)
if(Idefix_CUSTOM_EOS)
  set(Idefix_CUSTOM_EOS_FILE "eos_custom.hpp" CACHE FILEPATH "Custom equation of state source file")
endif()
set(Idefix_RECONSTRUCTION "Linear" CACHE STRING "Type of cell reconstruction scheme")
option(Idefix_HDF5 "Enable HDF5 I/O (requires HDF5 library)" OFF)
if(Idefix_MHD)
  option(Idefix_EVOLVE_VECTOR_POTENTIAL "Evolve the vector potential instead of the field (helps reducing div(B) in long runs)" OFF)
endif()
set_property(CACHE Idefix_RECONSTRUCTION PROPERTY STRINGS Constant Linear LimO3 Parabolic)
set(Idefix_PRECISION "Double" CACHE STRING "Precision of arithmetics")
set_property(CACHE Idefix_PRECISION PROPERTY STRINGS Double Single)

set(Idefix_LOOP_PATTERN "Default" CACHE STRING "Loop pattern for idefix_for")
set_property(CACHE Idefix_LOOP_PATTERN PROPERTY STRINGS Default SIMD Range MDRange TeamPolicy TeamPolicyInnerVector)


# load git revision tools
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)
include(ReplaceIdefixSource)
include(AddIdefixSource)
include(SetIdefixProperty)
include(SetRequiredBuildSettingsForGCC8)

#Idefix requires Cuda Lambdas (experimental)
if(Kokkos_ENABLE_CUDA)
  set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "Idefix requires lambdas on Cuda" FORCE)
endif()

# Add kokkos CMAKE files (required early since these set compiler options)
add_subdirectory(src/kokkos build/kokkos)

# Add Idefix CXX Flags
add_compile_options(${Idefix_CXX_FLAGS})

# Add filesystem libraries for GCC8
set_required_build_settings_for_GCC8()

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
  target_sources(idefix
    PUBLIC src/mpi.cpp
    PUBLIC src/mpi.hpp
  )
endif()

if(Idefix_HDF5)
  add_compile_definitions("WITH_HDF5")
  if(Idefix_MPI)
    set(HDF5_PREFER_PARALLEL TRUE)
  endif()
  target_sources(idefix
    PUBLIC src/output/xdmf.cpp
    PUBLIC src/output/xdmf.hpp
  )
  find_package(HDF5 REQUIRED)
  target_link_libraries(idefix "${HDF5_LIBRARIES}")
  target_include_directories(idefix PUBLIC "${HDF5_INCLUDE_DIRS}")
  message(STATUS "XDMF (hdf5+xmf) dumps enabled")
else()
  set(Idefix_HDF5 OFF)
endif()

if(Idefix_DEBUG)
  add_compile_definitions("DEBUG")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
endif()

if(Idefix_RUNTIME_CHECKS)
  add_compile_definitions("RUNTIME_CHECKS")
endif()

if(Idefix_HIGH_ORDER_FARGO)
  add_compile_definitions("HIGH_ORDER_FARGO")
endif()

if(Idefix_EVOLVE_VECTOR_POTENTIAL)
  add_compile_definitions("EVOLVE_VECTOR_POTENTIAL")
endif()
#update version.hpp if possible
git_describe(GIT_SHA1)
set(Idefix_VERSION ${Idefix_VERSION_MAJOR}.${Idefix_VERSION_MINOR}.${Idefix_VERSION_PATCH}-${GIT_SHA1})
file(WRITE ${CMAKE_SOURCE_DIR}/src/version.hpp "#define IDEFIX_GIT_COMMIT \"${GIT_SHA1}\"\n")
file(APPEND ${CMAKE_SOURCE_DIR}/src/version.hpp "#define IDEFIX_VERSION_MAJOR \"${Idefix_VERSION_MAJOR}\"\n")
file(APPEND ${CMAKE_SOURCE_DIR}/src/version.hpp "#define IDEFIX_VERSION_MINOR \"${Idefix_VERSION_MINOR}\"\n")
file(APPEND ${CMAKE_SOURCE_DIR}/src/version.hpp "#define IDEFIX_VERSION_PATCH \"${Idefix_VERSION_PATCH}\"\n")
file(APPEND ${CMAKE_SOURCE_DIR}/src/version.hpp "#define IDEFIX_VERSION \"${Idefix_VERSION}\"\n")


if(NOT ${Idefix_DEFS} STREQUAL "definitions.hpp")
  add_compile_definitions("DEFINITIONS_FILE=\"${Idefix_DEFS}\"")
endif()

if(Idefix_CUSTOM_EOS)
  add_compile_definitions("EOS_FILE=\"${Idefix_CUSTOM_EOS_FILE}\"")
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

#Loop type
if(${Idefix_LOOP_PATTERN} STREQUAL "SIMD")
  if(Kokkos_ENABLE_OPENMP)
    message(ERROR "SIMD loop pattern is incompatible with OPENMP")
  endif()
  if(Kokkos_ENABLE_CUDA)
    message(ERROR "SIMD loop pattern is incompatible with Cuda")
  endif()
  if(Kokkos_ENABLE_HIP)
    message(ERROR "SIMD loop pattern is incompatible with HIP")
  endif()
  add_compile_definitions("LOOP_PATTERN_SIMD")
elseif(${Idefix_LOOP_PATTERN} STREQUAL "Range")
  add_compile_definitions("LOOP_PATTERN_1DRANGE")
elseif(${Idefix_LOOP_PATTERN} STREQUAL "MDRange")
  add_compile_definitions("LOOP_PATTERN_MDRANGE")
elseif(${Idefix_LOOP_PATTERN} STREQUAL "TeamPolicy")
  add_compile_definitions("LOOP_PATTERN_TPX")
elseif(${Idefix_LOOP_PATTERN} STREQUAL "TeamPolicyInnerVector")
  add_compile_definitions("LOOP_PATTERN_TPTTRTVR")
elseif(NOT ${Idefix_LOOP_PATTERN} STREQUAL "Default")
  message(ERROR "Unknown loop Pattern")
endif()

# precision
if(${Idefix_PRECISION} STREQUAL "Single")
  add_compile_definitions("SINGLE_PRECISION")
endif()

target_include_directories(idefix PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )
target_include_directories(idefix PUBLIC
                           src/kokkos/core/src
                           src/dataBlock
                           src/dataBlock/planetarySystem
                           src/fluid
                           src/fluid/boundary
                           src/fluid/braginskii
                           src/fluid/constrainedTransport
                           src/fluid/eos
                           src/fluid/RiemannSolver
                           src/fluid/RiemannSolver/HDsolvers
                           src/fluid/RiemannSolver/MHDsolvers
                           src/fluid/RiemannSolver/Dustsolvers
                           src/fluid/tracer
                           src/output
                           src/rkl
                           src/gravity
                           src/utils
                           src/utils/iterativesolver
                           src
                           )

target_link_libraries(idefix Kokkos::kokkos)

message(STATUS "Idefix final configuration")
if(Idefix_EVOLVE_VECTOR_POTENTIAL)
  message(STATUS "    MHD:  ${Idefix_MHD} (Vector potential)")
else()
  message(STATUS "    MHD:  ${Idefix_MHD}")
endif()
message(STATUS "    MPI:  ${Idefix_MPI}")
message(STATUS "    HDF5: ${Idefix_HDF5}")
message(STATUS "    Reconstruction: ${Idefix_RECONSTRUCTION}")
message(STATUS "    Precision: ${Idefix_PRECISION}")
message(STATUS "    Version: ${Idefix_VERSION}")
message(STATUS "    Problem definitions: '${Idefix_DEFS}'")
if(Idefix_CUSTOM_EOS)
  message(STATUS "    EOS: Custom file '${Idefix_CUSTOM_EOS_FILE}'")
endif()
