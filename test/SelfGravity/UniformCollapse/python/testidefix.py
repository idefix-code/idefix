#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:29:41 2020

@author: glesur
"""

import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))
from pytools.vtk_io import readVTK
import numpy as np

# Arguments
gravCst = 1./(4.*np.pi) # grav cst 4piG
rho0 = 1e-3 # initial density

# Parameters (calculated from the arguments)
tff = np.sqrt(3.*np.pi/(32.*rho0*gravCst)) # Theoretical free-fall time to achieve collapse

# Collecting data
V=readVTK('../data.0005.vtk')
time = V.t[0]
rho = np.squeeze(V.data['RHO'])

# A function to detect density plateau, its size and position
def get_plateau(rho):
    plateau = [rho[0]]
    for dsty in rho[1:]:
        if np.abs(dsty-plateau[-1])>0.01*plateau[-1]:
            return np.mean(plateau)
        else:
            plateau.append(dsty)
    return "Error : the whole density distribution is approximately constant !"

# Isolating density plateau
plateau=get_plateau(rho)

# Calculating absolute error
r_th=time/tff # Theoretical ratio
eta = np.arccos((plateau/rho0)**(-1./6.))
r_num = 2./np.pi*(eta + 1./2.*np.sin(2.*eta)) # Numerical ratio
error = np.abs(r_num-r_th) # Absolute error

# Print the result of the test
print("Error=%e"%error)
if error<2.0e-3:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
