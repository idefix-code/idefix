######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/Radiation/novae_like_3D_no_damping_rad/CMakeLists.txt
#
# Last modified : 04/2025
#
# Copyright(C) by :
# - Nicolas Scepi <nicolas.scepi@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

add_idefix_source(analysis.cpp)
add_idefix_source(analysis.hpp)
enable_idefix_property(Idefix_MHD)
