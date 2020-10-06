#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:29:41 2020

@author: glesur
"""

import idefixTools as idfx
import numpy as np
import sys
import argparse
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

parser = argparse.ArgumentParser()
parser.add_argument("-noplot",
                    default=False,
                    help="disable plotting",
                    action="store_true")


args=parser.parse_args()

V=idfx.readVTKCart('../data.0001.vtk')
U=idfx.readVTKCart('data.ref.vtk')

solinterp=interp1d(U.x,U.data['PRS'][:,0,0])

if(not args.noplot):
    plt.figure(1)
    plt.plot(V.x,V.data['RHO'][:,0,0],'+',markersize=2)
    plt.plot(U.x,U.data['RHO'][:,0,0])
    plt.title('Density')

    plt.figure(2)
    #plt.plot(x,u)
    plt.plot(V.x,V.data['VX1'][:,0,0],'+',markersize=2)
    plt.plot(U.x,U.data['VX1'][:,0,0])
    plt.title('X Velocity')

    plt.figure(3)
    #plt.plot(x,u)
    plt.plot(V.x,V.data['VX2'][:,0,0],'+',markersize=2)
    plt.plot(U.x,U.data['VX2'][:,0,0])
    plt.title('Y Velocity')

    plt.figure(4)
    #plt.plot(x,u)
    plt.plot(V.x,V.data['BX1'][:,0,0],'+',markersize=2)
    plt.plot(U.x,U.data['BX1'][:,0,0])
    plt.title('X field')

    plt.figure(5)
    #plt.plot(x,u)
    plt.plot(V.x,V.data['BX2'][:,0,0],'+',markersize=2)
    plt.plot(U.x,U.data['BX2'][:,0,0])
    plt.title('Y field')


    plt.figure(6)
    #plt.plot(x,p)
    plt.plot(V.x,V.data['PRS'][:,0,0],'+',markersize=2)
    #plt.plot(V.x,solinterp(V.x))
    plt.plot(U.x,U.data['PRS'][:,0,0])
    plt.title('Pressure')

    plt.ioff()
    plt.show()

error=np.mean(np.fabs(V.data['PRS'][:,0,0]-solinterp(V.x))**2/solinterp(V.x))

print("Error=%e"%error)
if error<1e-2:
    print("SUCCESS!")
else:
    print("FAILURE!")
