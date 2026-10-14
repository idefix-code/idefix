######################################################################################
# Idefix MHD astrophysical code
#
# Source file cmake/HandleDeps.cmake
#
# Last modified : 10/2026
#
# Copyright(C) by :
# - Sébastien Valat <sebastien.valat@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

function(GoogleTestDownloadAndLink)
  # CAUTION: if you change this verison, you should also check if you need to adapt
  # the fix in gtest-idfx and re-import the definition of GTEST_TEST_().
  set(gtest_version "1.18.0")

  # it not, download it internally to the project
  cmake_policy(SET CMP0135 NEW)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/releases/download/v${gtest_version}/googletest-${gtest_version}.tar.gz
  )
  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
endfunction()
