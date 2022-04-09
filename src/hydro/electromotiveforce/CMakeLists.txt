target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/calcCornerEmf.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/calcNonidealEMF.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/calcRiemannEmf.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/electroMotiveForce.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/electroMotiveForce.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/EMFexchange.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/enforceEMFBoundary.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/evolveMagField.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/evolveVectorPotential.cpp
  )
