#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 15:52:44 2021

@author: lesurg
"""
import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))
from pytools.idfx_io import readIdfxFile
import matplotlib.pyplot as plt

rep="../"

idfx=readIdfxFile(rep+"analysis.99.0.idfx")
d=idfx.data
Bx1s=d['Vs'][0,:-1,:-1,:]
Bx2s=d['Vs'][1,:-1,:,:-1]
Bx3s=d['Vs'][2,:,:-1,:-1]

nxhalf=10
Bx1slice=Bx1s[:,:,nxhalf]
Bx2slice=Bx2s[:,:,nxhalf]
Bx3slice=Bx3s[:,:,nxhalf]


plt.matshow(Bx1slice,cmap='RdBu_r')
plt.colorbar()
plt.title("Bx1")

plt.matshow(Bx2slice,cmap='RdBu_r')
plt.colorbar()
plt.title("Bx2")
plt.matshow(Bx3slice,cmap='RdBu_r')
plt.colorbar()
plt.title("Bx3")


# inner boundary
Bx1slice=Bx1s[:,:5,nxhalf]
Bx2slice=Bx2s[:,:5,nxhalf]
Bx3slice=Bx3s[:,:5,nxhalf]

Vx1slice=d['Vc'][1,:,:5,nxhalf]
Vx2slice=d['Vc'][2,:,:5,nxhalf]
Vx3slice=d['Vc'][3,:,:5,nxhalf]


plt.matshow(Bx1slice,cmap='RdBu_r')
plt.colorbar()
plt.title("Bx1")
plt.matshow(Bx2slice,cmap='RdBu_r')
plt.colorbar()
plt.title("Bx2")
plt.matshow(Bx3slice,cmap='RdBu_r')
plt.colorbar()
plt.title("Bx3")

plt.matshow(Vx1slice,cmap='RdBu_r')
plt.colorbar()
plt.title("Vx1")
plt.matshow(Vx2slice,cmap='RdBu_r')
plt.colorbar()
plt.title("Vx2")
plt.matshow(Vx3slice,cmap='RdBu_r')
plt.colorbar()
plt.title("Vx3")

plt.show()
