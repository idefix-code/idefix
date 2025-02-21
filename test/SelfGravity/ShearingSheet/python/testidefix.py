#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:42:19 2025

@author: vdmba
"""
import os
import sys
import numpy as np
import inifix
sys.path.append(os.getenv("IDEFIX_DIR"))
from pytools.idfx_io import readIdfxFile
from pytools.vtk_io import readVTK
import matplotlib.pyplot as plt
import argparse

# check that the density the potential minima is sheared
# at the correct rate in the boundaries


parser = argparse.ArgumentParser()
parser.add_argument("-noplot",
                    default=True,
                    help="disable plotting",
                    action="store_true")


args, unknown=parser.parse_known_args()


Vbeg = readVTK("../data.0000.vtk")
Vend = readVTK("../data.0002.vtk")

Iend=readIdfxFile('../test.2.0.idfx') # Last Idfx
Iref=readIdfxFile('../test.0.0.idfx') # Ref vtk (initial state)

Ly = Vbeg.yl[-1]-Vbeg.yl[0]
xLeft = Vbeg.xl[0]
xRight = Vbeg.xl[-1]

deltaT = Vend.t[0]-Vbeg.t[0]

conf = inifix.load("../idefix.ini")

S = conf["Hydro"]["shearingBox"]

RHO = 0
nghost = 2
dy = Ly/Vbeg.ny


indexMaxInitLeftRHO = np.argmax(Iref.data["Vc"][RHO,0,nghost:-nghost,0])
indexMaxInitRightRHO = np.argmax(Iref.data["Vc"][RHO,0,nghost:-nghost,-1])

posMaxInitLeftRHO = Vbeg.y[indexMaxInitLeftRHO]
posMaxInitRightRHO = Vbeg.y[indexMaxInitRightRHO]

#breakpoint()
if args.noplot:
    fig, ax = plt.subplots()
    ax.plot(Vbeg.y,Iend.data["Pot"][RHO,0,nghost:-nghost,0]/Iend.data["Pot"][RHO,0,nghost:-nghost,0].max())
    ax.plot(Vbeg.y,Iend.data["Pot"][RHO,0,nghost:-nghost,-1]/Iend.data["Pot"][RHO,0,nghost:-nghost,-1].max())
    plt.show()

# potential is not computed at 0th step: using the position of the initial density maximum
posMaxInitLeftPot = posMaxInitLeftRHO
posMaxInitRightPot = posMaxInitRightRHO

indexMaxEndLeftPot = np.argmin(Iend.data["Pot"][RHO,0,nghost:-nghost,0])
indexMaxEndRightPot = np.argmin(Iend.data["Pot"][RHO,0,nghost:-nghost,-1])

posMaxEndLeftPot = Vbeg.y[indexMaxEndLeftPot]
posMaxEndRightPot = Vbeg.y[indexMaxEndRightPot]

deltaPosLeftPot = posMaxEndLeftPot - posMaxInitLeftPot
deltaPosRightPot = posMaxEndRightPot - posMaxInitRightPot

errorLeftPot = abs(deltaPosLeftPot - S*xLeft*deltaT)
errorRightPot = abs(deltaPosRightPot - S*xRight*deltaT)

errorPot = (errorLeftPot**2 + errorRightPot**2)**.5


err=abs(1-errorPot/dy)
print("Error=",err)

if(err<0.05):
  print("SUCCESS")
  sys.exit(0)
else:
  print("FAILED")
  sys.exit(1)
