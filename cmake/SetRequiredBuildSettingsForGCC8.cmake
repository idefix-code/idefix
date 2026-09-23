######################################################################################
# Idefix MHD astrophysical code
#
# Source file cmake/SetRequiredBuildSettingsForGCC8.cmake
#
# Last modified : 11/2023
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2023)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

# The std::filesystem functions are part of the C++17 standard, but GCC8
# does not store them in libstdc++ but in libstdc++-fs. This function
# ensure that we link to this additional library when using GCC8
# from https://discourse.cmake.org/t/correct-way-to-link-std-filesystem-with-gcc-8/4121
function( set_required_build_settings_for_GCC8 )
    # Always link with libstdc++fs.a when using GCC 8.
    # Note: This command makes sure that this option comes pretty late on the cmdline.
    link_libraries( "$<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:-lstdc++fs>" )
endfunction()
