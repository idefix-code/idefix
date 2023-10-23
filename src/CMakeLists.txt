
add_subdirectory(fluid)
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
# NOTE: mpi file are added by the root CMakeLists along with dependencies to the MPI library


#add_subdirectory("${PROJECT_BINARY_DIR}" "${PROJECT_BINARY_DIR}/local")
