target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/dump.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/dump.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/output.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/output.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/vtk.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/vtk.hpp
  )
if(Idefix_HDF5)
target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/xdmf.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/xdmf.hpp
  )
endif()
