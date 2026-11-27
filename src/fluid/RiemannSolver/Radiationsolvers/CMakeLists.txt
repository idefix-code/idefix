######################################################################################
# Idefix MHD astrophysical code
#
# Source file src/fluid/RiemannSolver/Radiationsolvers/CMakeLists.txt
#
# Last modified : 02/2026
#
# Copyright(C) by :
# - Nicolas Scepi <nicolas.scepi@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

target_sources(idefix
PUBLIC ${CMAKE_CURRENT_LIST_DIR}/hllRad.hpp
PUBLIC ${CMAKE_CURRENT_LIST_DIR}/hllcRad.hpp
PUBLIC ${CMAKE_CURRENT_LIST_DIR}/lfrRad.hpp
PUBLIC ${CMAKE_CURRENT_LIST_DIR}/speedRad.hpp
PUBLIC ${CMAKE_CURRENT_LIST_DIR}/lim_fluxRad.hpp
)
