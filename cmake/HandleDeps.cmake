######################################################################################
# Idefix MHD astrophysical code
#
# Source file cmake/AddIdefixSource.cmake
#
# Last modified : 12/2021
#
# Copyright(C) by :
# - Sébastien Valat <sebastien.valat@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

function(FindSystemGoogleTestOrDownload)
  # try to use the system one
  find_package(GTest QUIET)

  # it not, download it internally to the project
  if (NOT GTest_FOUND)
    cmake_policy(SET CMP0135 NEW)
    include(FetchContent)
    FetchContent_Declare(
      googletest
      URL https://github.com/google/googletest/releases/download/v1.18.0/googletest-1.18.0.tar.gz
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
  endif()
endfunction()
