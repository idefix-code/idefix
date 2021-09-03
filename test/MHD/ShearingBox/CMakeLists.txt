# add specific source file to idefix build tree
target_sources(idefix
    PUBLIC ${CMAKE_CURRENT_LIST_DIR}/analysis.cpp
    PUBLIC ${CMAKE_CURRENT_LIST_DIR}/analysis.hpp
)
