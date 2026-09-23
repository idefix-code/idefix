######################################################################################
# Idefix MHD astrophysical code
#
# Source file cmake/AddIdefixSource.cmake
#
# Last modified : 12/2021
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2021)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

function(add_idefix_source _new_source)
  message(STATUS "    Adding problem-specific source file ${_new_source}")
  target_sources(idefix PUBLIC ${_new_source})
endfunction()
