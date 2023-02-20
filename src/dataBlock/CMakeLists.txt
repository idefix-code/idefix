add_subdirectory(planetarySystem)

target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/coarsen.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/dataBlock.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/dataBlock.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/dataBlockHost.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/dataBlockHost.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/dumpToFile.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/evolveStage.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/fargo.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/fargo.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/makeGeometry.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/stateContainer.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/stateContainer.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/validation.cpp

  )
