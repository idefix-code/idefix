######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/MHD/ShearingBox/CMakeLists.txt
#
# Last modified : 12/2021
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2021)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

# add specific source file to idefix build tree
add_idefix_source(analysis.cpp)
add_idefix_source(analysis.hpp)
enable_idefix_property(Idefix_MHD)
