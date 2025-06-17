#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:42:19 2025

@author: vdbma
"""
import os
import sys
import numpy as np
import inifix
sys.path.append(os.getenv("IDEFIX_DIR"))
from pytools.vtk_io import readVTK
import matplotlib.pyplot as plt
import argparse

# check that the potential minima is sheared
# at the correct rate in the boundaries


parser = argparse.ArgumentParser()
parser.add_argument("-noplot",
                    default=True,
                    help="disable plotting",
                    action="store_true")


args, unknown=parser.parse_known_args()


Vbeg = readVTK("../data.0000.vtk")
Vend = readVTK("../data.0002.vtk")

Ly = Vbeg.yl[-1]-Vbeg.yl[0]
xLeft = Vbeg.x[0]
xRight = Vbeg.x[-1]

deltaT = Vend.t[0]-Vbeg.t[0]

conf = inifix.load("../idefix.ini")

S = conf["Hydro"]["shearingBox"]

dy = Ly/Vbeg.ny


indexMaxInitLeftRHO = np.argmax(Vbeg.data["RHO"][0,:,0])
indexMaxInitRightRHO = np.argmax(Vbeg.data["RHO"][-1,:,0])

posMaxInitLeftRHO = Vbeg.y[indexMaxInitLeftRHO]
posMaxInitRightRHO = Vbeg.y[indexMaxInitRightRHO]

#breakpoint()
if not args.noplot:
    fig, ax = plt.subplots()
    ax.plot(Vend.y,Vend.data["phiP"][0,:,0]/Vend.data["phiP"][0,:,:].max(),".")
    ax.plot(Vend.y,Vend.data["phiP"][-1,:,0]/Vend.data["phiP"][-1,:,0].max(),".")
    plt.show()

# potential is not computed at 0th step: using the position of the initial density maximum
posMaxInitLeftPot = posMaxInitLeftRHO
posMaxInitRightPot = posMaxInitRightRHO

indexMaxEndLeftPot = np.argmin(Vend.data["phiP"][0,:,0])
indexMaxEndRightPot = np.argmin(Vend.data["phiP"][-1,:,0])

posMaxEndLeftPot = Vbeg.y[indexMaxEndLeftPot]
posMaxEndRightPot = Vbeg.y[indexMaxEndRightPot]

deltaPosLeftPot = posMaxEndLeftPot - posMaxInitLeftPot
deltaPosRightPot = posMaxEndRightPot - posMaxInitRightPot

errorLeftPot = abs(deltaPosLeftPot - S*xLeft*deltaT)
errorRightPot = abs(deltaPosRightPot - S*xRight*deltaT)

errorPot = (errorLeftPot**2 + errorRightPot**2)**.5



err=abs(errorPot/dy)
print("Error=",err)

if(err<1): # total cummulated shift is less than on cell
  print("SUCCESS")
  sys.exit(0)
else:
  print("FAILED")
  sys.exit(1)
