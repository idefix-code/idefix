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
from scipy.fftpack import fft, fftfreq

# Arguments
gravCst = np.pi # grav cst G=pi
gamma = 5/3
rho_bg = 3.0 # background density

# Parameters (calculated from the arguments)
prs_bg = 1.0/gamma
cs02 = gamma*prs_bg/rho_bg # squared bg speed velocity
lbdJ = np.sqrt(cs02*np.pi/(gravCst*rho_bg)) # Jean's wavelength
kJ = 2*np.pi/lbdJ # Jean's wavevector

# Collect vtk files (we ignore those before vtkRef)
Vmax=readVTK('../data.0004.vtk') # Best time vtk (for gr fitting)
Vref=readVTK('../data.0003.vtk') # Ref vtk (after unwanted modes are quenched)

# Collect data
x1 = np.squeeze(Vref.x)
tref, tmax = Vref.t[0], Vmax.t[0]
uref, umax = np.squeeze(Vref.data['VX1']), np.squeeze(Vmax.data['VX1'])

# Compute fft
N1 = x1.size
dx1 = x1[1]-x1[0]
k1 = 2*np.pi*fftfreq(N1, dx1)
furef, fumax = np.abs(fft(uref)), np.abs(fft(umax))

# Compute the numerical/theoretical growth rate of the second non-zero unstable mode
mode = 2
gr_num = np.log(fumax/furef)[mode]/(tmax-tref)
gr_th = np.sqrt(np.abs(k1[mode]**2*cs02-4.*np.pi*gravCst*rho_bg))
error = np.abs((gr_num-gr_th)/gr_th)

# Print the result of the test
print("Error=%e"%error)
if error<2.3e-3:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
