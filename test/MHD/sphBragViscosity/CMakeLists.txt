######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/MHD/sphBragViscosity/CMakeLists.txt
#
# Last modified : 10/2023
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG / CNRS - 2023)
# - and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

enable_idefix_property(Idefix_MHD)
## The current test is not physics, but merely a non-regression check
## The physical test needs to have the divergence of the velocity field artifiacially set to zero
#replace_idefix_source(fluid/braginskii/bragViscosity.cpp noDivBragViscosity.cpp)
