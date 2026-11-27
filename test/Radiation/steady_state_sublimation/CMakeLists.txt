######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/Radiation/steady_state_sublimation/CMakeLists.txt
#
# Last modified : 09/2025
#
# Copyright(C) by :
# - Nicolas Scepi <nicolas.scepi@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

replace_idefix_source(dataBlock/evolveStage.cpp evolveStage.cpp)
replace_idefix_source(fluid/radsource.cpp radsource.cpp)
