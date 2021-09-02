
add_subdirectory(hydro)
add_subdirectory(dataBlock)
add_subdirectory(output)
add_subdirectory(rkl)

set(idefix_source_file ${idefix_source_file}
  ${CMAKE_CURRENT_LIST_DIR}/arrays.hpp
  ${CMAKE_CURRENT_LIST_DIR}/error.cpp
  ${CMAKE_CURRENT_LIST_DIR}/error.hpp
  ${CMAKE_CURRENT_LIST_DIR}/global.cpp
  ${CMAKE_CURRENT_LIST_DIR}/global.hpp
  ${CMAKE_CURRENT_LIST_DIR}/grid.cpp
  ${CMAKE_CURRENT_LIST_DIR}/grid.hpp
  ${CMAKE_CURRENT_LIST_DIR}/gridHost.hpp
  ${CMAKE_CURRENT_LIST_DIR}/gridHost.cpp
  ${CMAKE_CURRENT_LIST_DIR}/idefix.hpp
  ${CMAKE_CURRENT_LIST_DIR}/input.cpp
  ${CMAKE_CURRENT_LIST_DIR}/input.hpp
  ${CMAKE_CURRENT_LIST_DIR}/loop.hpp
  ${CMAKE_CURRENT_LIST_DIR}/macros.hpp
  ${CMAKE_CURRENT_LIST_DIR}/main.cpp
  ${CMAKE_CURRENT_LIST_DIR}/profiler.cpp
  ${CMAKE_CURRENT_LIST_DIR}/profiler.hpp
  ${CMAKE_CURRENT_LIST_DIR}/real_types.hpp
  ${CMAKE_CURRENT_LIST_DIR}/setup.hpp
  ${CMAKE_CURRENT_LIST_DIR}/timeIntegrator.hpp
  ${CMAKE_CURRENT_LIST_DIR}/timeIntegrator.cpp
  PARENT_SCOPE
  )


#add_subdirectory("${PROJECT_BINARY_DIR}" "${PROJECT_BINARY_DIR}/local")
