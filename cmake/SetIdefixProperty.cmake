######################################################################################
# Idefix MHD astrophysical code
#
# Source file cmake/SetIdefixProperty.cmake
#
# Last modified : 12/2021
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2021)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

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
