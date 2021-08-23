
add_subdirectory(kokkos)
add_subdirectory(hydro)
add_subdirectory(dataBlock)
add_subdirectory(output)
add_subdirectory(rkl)

set(idefix_source_file ${idefix_source_file}
  src/arrays.hpp
  src/error.cpp
  src/error.hpp
  src/global.cpp
  src/global.hpp
  src/grid.cpp
  src/grid.hpp
  src/gridHost.hpp
  src/gridHost.cpp
  src/idefix.hpp
  src/input.cpp
  src/input.hpp
  src/loop.hpp
  src/macros.hpp
  src/main.cpp
  src/profiler.cpp
  src/profiler.hpp
  src/real_types.hpp
  src/setup.hpp
  src/timeIntegrator.hpp
  src/timeIntegrator.cpp
  PARENT_SCOPE
  )


#add_subdirectory("${PROJECT_BINARY_DIR}" "${PROJECT_BINARY_DIR}/local")
