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
import argparse
import numpy as np
import matplotlib.pyplot as plt
import inifix

parser = argparse.ArgumentParser()
parser.add_argument("-noplot",
                    default=False,
                    help="disable plotting",
                    action="store_true")


args, unknown=parser.parse_known_args()

# find which ini file was used
with open('../idefix.0.log','r') as fp:
  lines = fp.readlines()
  for line in lines:
    if line.find(".ini") != -1:
      s=line.split()
      inputFile=s[5][:-1]

input = inifix.load('../'+inputFile)
componentDir = input["Setup"]["direction"][0]
spatialDir = input["Setup"]["direction"][1]

# Convert componentDir
componentName=["BX1","BX2","BX3"]


V=readVTK('../data.0001.vtk', geometry='cartesian')

# left_state and right_state set p, rho and u
# geometry sets left boundary on 0., right boundary on 1 and initial
# position of the shock xi on 0.5
# t is the time evolution for which positions and states in tube should be calculated
# gamma denotes specific heat
# note that gamma and npts are default parameters (1.4 and 500) in solve function

x=V.x
if spatialDir==0:
  qside=V.data[componentName[componentDir]][:,4,0]
  qcenter=V.data[componentName[componentDir]][:,32,0]
  x=V.x
if spatialDir==1:
  qside=V.data[componentName[componentDir]][4,:,0]
  qcenter=V.data[componentName[componentDir]][32,:,0]
  x=V.y
if spatialDir==2:
  qside=V.data[componentName[componentDir]][4,0,:]
  qcenter=V.data[componentName[componentDir]][32,0,:]
  x=V.z

#initial condition in the code
q0=0.*x+0.1*(np.abs(x)<0.1)

# compute solution to diffusion equation using Fourier transforms
qf0=np.fft.fft(q0)
k2=(2.0*np.pi*q0.size*np.fft.fftfreq(q0.size))**2
qf=qf0*np.exp(-k2*V.t*0.05)
qtheory=np.real(np.fft.ifft(qf))

if(not args.noplot):
    plt.figure(1)
    plt.plot(x,qside,'+',markersize=6,label='side')
    plt.plot(x,qcenter,'+',markersize=6,label='center')

    #plt.plot(x,q0)
    plt.plot(x,qtheory,label='theory')
    plt.legend()
    plt.title(componentName[componentDir])



    plt.ioff()
    plt.show()

errorCenter=np.mean((qcenter-qtheory)**2)
errorSide=np.mean((qside-qtheory)**2)
print("Error Side=%e, Error Center=%e"%(errorSide,errorCenter))
assert(errorCenter<2e-5)
assert(errorSide<2e-8)

print("SUCCESS!")
