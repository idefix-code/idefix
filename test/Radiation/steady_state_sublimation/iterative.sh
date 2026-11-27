#!/bin/bash
######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/Radiation/steady_state_sublimation/iterative.sh
#
# Last modified : 09/2025
#
# Copyright(C) by :
# - Nicolas Scepi <nicolas.scepi@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

cp 'dump.hr002.dmp' 'dump.ini.dmp'

for i in {0..15}
do
  mpirun --mca io ^ompio -np 2 ./idefix -dec 2 1
  cp 'dump.0001.dmp' 'dump.ini.dmp'
  cp 'dump.0001.dmp' 'dump.cycle'$i'.dmp'
  mkdir 'cycle'$i
  mv data.*.vtk 'cycle'$i/
done
