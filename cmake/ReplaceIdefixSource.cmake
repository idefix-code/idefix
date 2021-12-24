function(replace_idefix_source _old_source _new_source)
  get_target_property(fullsource idefix SOURCES)
  list(FILTER fullsource INCLUDE REGEX ${_old_source})

  # Check that only one file matches the pattern
  list(LENGTH fullsource listsize)
  if(listsize GREATER 1)
    message(SEND_ERROR "several source files match your search pattern ${_old_source}")
  endif()
  if(listsize LESS 1)
    message(SEND_ERROR "no source file match your search pattern ${_old_source}")
  endif()
  message(STATUS "    Replacing: ${fullsource} by ${_new_source}")

  # remove the old source file from the target source list
  get_target_property(mylist idefix SOURCES)
  list(FILTER mylist EXCLUDE REGEX ${_old_source})
  # replace the source with our new source list
  set_property(TARGET idefix PROPERTY SOURCES ${mylist})
  # add the new source file
  target_sources(idefix PUBLIC ${_new_source})
endfunction()
