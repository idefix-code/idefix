function(set_idefix_property _property_name _property_value)
  message(STATUS "    Explicitely setting ${_property_name} = ${_property_value}")
  set(${_property_name} ${_property_value} CACHE STRING "Overridden by custom setup in CMakeLists.txt" FORCE)
endfunction()

function(enable_idefix_property _property_name)
  message(STATUS "    Explicitely setting ${_property_name} = ON")
  set(${_property_name} ON CACHE STRING "Overridden by custom setup in CMakeLists.txt" FORCE)
endfunction()

function(disable_idefix_property _property_name)
  message(STATUS "    Explicitely setting ${_property_name} = OFF")
  set(${_property_name} OFF CACHE STRING "Overridden by custom setup in CMakeLists.txt" FORCE)
endfunction()
