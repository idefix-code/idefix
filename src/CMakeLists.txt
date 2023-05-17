
add_subdirectory(hydro)
add_subdirectory(dataBlock)
add_subdirectory(output)
add_subdirectory(rkl)
add_subdirectory(gravity)
add_subdirectory(utils)

target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/arrays.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/error.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/error.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/global.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/global.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/grid.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/grid.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/gridHost.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/gridHost.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/idefix.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/input.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/input.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/loop.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/macros.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/main.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/profiler.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/profiler.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/real_types.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/reduce.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/setup.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/setup.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/timeIntegrator.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/timeIntegrator.cpp
  )

if(Idefix_MPI)
  target_sources(idefix
    PUBLIC ${CMAKE_CURRENT_LIST_DIR}/mpi.cpp
    PUBLIC ${CMAKE_CURRENT_LIST_DIR}/mpi.hpp
    )
endif()

if(Idefix_HDF5)
  if(Idefix_MPI)
    set(HDF5_PREFER_PARALLEL TRUE)
  endif()
  find_package(HDF5 REQUIRED)
  target_link_libraries(idefix "${HDF5_LIBRARIES}")
  target_include_directories(idefix PUBLIC "${HDF5_INCLUDE_DIRS}")
  message(STATUS "HIGHFIVE Found: ${CMAKE_CURRENT_LIST_DIR}/HighFive/include")
  message(STATUS "XDMF (hdf5+xmf) dumps enabled")
  target_include_directories(idefix PUBLIC "${CMAKE_CURRENT_LIST_DIR}/HighFive/include")
endif()

#add_subdirectory("${PROJECT_BINARY_DIR}" "${PROJECT_BINARY_DIR}/local")
