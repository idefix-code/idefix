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

V=readVTK('../data.0001.vtk')
R=readVTK('1Dsolution/data.0001.vtk')



# Finally, let's plot solutions
p = R.data['PRS'][:,0,0]
x= R.r

[x3D, y3D, z3D]= np.meshgrid(V.x,V.y,V.z,indexing='ij')
r3D=np.sqrt(x3D**2+y3D**2+z3D**2)
p3D=V.data['PRS']


solinterp=interp1d(x,p)


if(not args.noplot):
    plt.figure(1)
    plt.plot(x,p)
    plt.plot(r3D.flatten(),p3D.flatten(),'+',markersize=2)
    plt.title('Pressure')

    plt.ioff()
    plt.show()

error=np.mean(np.fabs(V.data['PRS'][:,0,0]-solinterp(V.x)))
print("Error=%e"%error)
if error<2e-3:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
