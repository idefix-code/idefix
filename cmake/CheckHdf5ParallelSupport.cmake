include(CheckCSourceCompiles)

# Check whether the HDF5 library found by find_package(HDF5) supports parallel
# (MPI-IO) access, and store the result in IDEFIX_HDF5_IS_PARALLEL in the
# caller's scope.
#
# Usage: CheckHdf5ParallelSupport(<hdf5_link_items>)
# where <hdf5_link_items> is the list of HDF5 targets/libraries to link against.
#
# Some CMake/HDF5 combinations do not reliably set HDF5_IS_PARALLEL (e.g. when
# HDF5 is found via its CMake config package rather than the FindHDF5 module).
# Combine metadata checks with a compile+link probe for the MPI-IO symbols to
# reliably detect parallel HDF5.
function(CheckHdf5ParallelSupport hdf5_link_items)
  set(_idefix_hdf5_is_parallel FALSE)
  if(HDF5_IS_PARALLEL OR HDF5_C_IS_PARALLEL)
    set(_idefix_hdf5_is_parallel TRUE)
  endif()

  if(NOT _idefix_hdf5_is_parallel AND DEFINED HDF5_C_COMPILER_EXECUTABLE)
    execute_process(
      COMMAND "${HDF5_C_COMPILER_EXECUTABLE}" -showconfig
      OUTPUT_VARIABLE _idefix_hdf5_showconfig
      ERROR_QUIET
    )
    if(_idefix_hdf5_showconfig MATCHES "Parallel HDF5:[ \t]*yes")
      set(_idefix_hdf5_is_parallel TRUE)
    endif()
  endif()

  set(CMAKE_REQUIRED_INCLUDES ${HDF5_INCLUDE_DIRS} ${MPI_C_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${hdf5_link_items} MPI::MPI_C)
  check_c_source_compiles(
    "#include <mpi.h>
     #include <hdf5.h>
     int main(void) {
       hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
       H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL);
       H5Pclose(plist);
       return 0;
     }"
    IDEFIX_HDF5_HAS_MPI_IO
  )

  if(IDEFIX_HDF5_HAS_MPI_IO)
    set(_idefix_hdf5_is_parallel TRUE)
  endif()

  set(IDEFIX_HDF5_IS_PARALLEL ${_idefix_hdf5_is_parallel} PARENT_SCOPE)
endfunction()
