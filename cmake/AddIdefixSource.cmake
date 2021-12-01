function(add_idefix_source _new_source)
  message(STATUS "    Adding problem-specific source file ${_new_source}")
  target_sources(idefix PUBLIC ${_new_source})
endfunction()
