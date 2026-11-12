target_sources(idefix-core
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/fft.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/fft.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/transpose.hpp
  )

if (Idefix_GTEST)
  target_sources(unittests
    PRIVATE ${CMAKE_CURRENT_LIST_DIR}/tests/testFFT.cpp
  )
endif()
