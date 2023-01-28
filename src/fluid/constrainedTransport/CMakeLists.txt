target_sources(idefix
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/calcCornerEmf.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/calcNonidealEMF.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/calcRiemannEmf.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/constrainedTransport.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/EMFexchange.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/enforceEMFBoundary.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/evolveMagField.hpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/evolveVectorPotential.hpp
  )
