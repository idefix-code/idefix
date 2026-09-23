######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/MHD/MTI/CMakeLists.txt
#
# Last modified : 10/2023
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2023)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

set_idefix_property(Idefix_RECONSTRUCTION LimO3)
enable_idefix_property(Idefix_MHD)
add_idefix_source(analysis.hpp)
add_idefix_source(analysis.cpp)
