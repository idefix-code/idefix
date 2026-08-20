target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/gravity.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/gravity.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/laplacian.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/laplacian.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/selfGravity.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/selfGravity.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/selfGravityIterative.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/selfGravityIterative.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/selfGravityFFT.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/selfGravityFFT.cpp
  )
