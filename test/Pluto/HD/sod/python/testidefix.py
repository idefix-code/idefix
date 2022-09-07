#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:29:41 2020

@author: glesur
"""

import idefixTools as idfx
import sod
import sys
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

V=idfx.readVTKCart('../data.0002.vtk')
gamma = 1.4
npts = 5000

# left_state and right_state set p, rho and u
# geometry sets left boundary on 0., right boundary on 1 and initial
# position of the shock xi on 0.5
# t is the time evolution for which positions and states in tube should be calculated
# gamma denotes specific heat
# note that gamma and npts are default parameters (1.4 and 500) in solve function
positions, regions, values = sod.solve(left_state=(1, 1, 0), right_state=(0.1, 0.125, 0.),
                                       geometry=(0., 1., 0.5), t=0.2, gamma=gamma, npts=npts)


# Finally, let's plot solutions
p = values['p']
rho = values['rho']
u = values['u']
x= values['x']


solinterp=interp1d(x,p)


if(not args.noplot):
    plt.figure(1)
    plt.plot(x,rho)
    plt.plot(V.x,V.data['rho'][:,0,0],'+',markersize=2)
    plt.title('Density')

    plt.figure(2)
    plt.plot(x,u)
    plt.plot(V.x,V.data['vx1'][:,0,0],'+',markersize=2)
    plt.title('Velocity')

    plt.figure(3)
    plt.plot(x,p)
    plt.plot(V.x,V.data['prs'][:,0,0],'+',markersize=2)
    plt.title('Pressure')

    plt.ioff()
    plt.show()

error=np.mean(np.fabs(V.data['prs'][:,0,0]-solinterp(V.x)))
print("Error=%e"%error)
if error<1e-2:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
