add_subdirectory(iterativesolver)

target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/bigEndian.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/dumpImage.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/dumpImage.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/lookupTable.hpp
  )
