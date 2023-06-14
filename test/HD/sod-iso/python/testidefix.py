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
from pytools import sod
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

parser = argparse.ArgumentParser()
parser.add_argument("-noplot",
                    default=False,
                    help="disable plotting",
                    action="store_true")


args, unknown=parser.parse_known_args()

V=readVTK('../data.0002.vtk', geometry='cartesian')
gamma = 1.00000000001
npts = 5000

# left_state and right_state set p, rho and u
# geometry sets left boundary on 0., right boundary on 1 and initial
# position of the shock xi on 0.5
# t is the time evolution for which positions and states in tube should be calculated
# gamma denotes specific heat
# note that gamma and npts are default parameters (1.4 and 500) in solve function
positions, regions, values = sod.solve(left_state=(1, 1, 0), right_state=(0.125, 0.125, 0.),
                                       geometry=(0., 1., 0.5), t=0.2, gamma=gamma, npts=npts)


# Finally, let's plot solutions
p = values['p']
rho = values['rho']
u = values['u']
x= values['x']


solinterp=interp1d(x,u)

error=np.mean(np.fabs((V.data['VX1'][:,0,0]-solinterp(V.x))))




if(not args.noplot):
    plt.figure(1)
    plt.plot(x,rho)
    plt.plot(V.x,V.data['RHO'][:,0,0],'+',markersize=2)
    plt.title('Density')

    plt.figure(2)
    plt.plot(x,u)
    plt.plot(V.x,V.data['VX1'][:,0,0],'+',markersize=2)
    plt.title('Velocity')

    plt.ioff()
    plt.show()

print("Error=%e"%error)
if error<5e-3:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
