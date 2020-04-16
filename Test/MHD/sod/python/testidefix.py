#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:29:41 2020

@author: glesur
"""

import idefixTools as idfx
import numpy as np
import matplotlib.pyplot as plt

V=idfx.readVTKCart('../data.0001.vtk')
U=idfx.readVTKCart('data.ref.vtk')





plt.figure(1)
plt.plot(V.x,V.data['Vc0'][:,0,0],'+',markersize=2)
plt.plot(U.x,U.data['Vc0'][:,0,0])
plt.title('Density')

plt.figure(2)
#plt.plot(x,u)
plt.plot(V.x,V.data['Vc1'][:,0,0],'+',markersize=2)
plt.plot(U.x,U.data['Vc1'][:,0,0])
plt.title('X Velocity')

plt.figure(3)
#plt.plot(x,u)
plt.plot(V.x,V.data['Vc2'][:,0,0],'+',markersize=2)
plt.plot(U.x,U.data['Vc2'][:,0,0])
plt.title('Y Velocity')

plt.figure(4)
#plt.plot(x,u)
plt.plot(V.x,V.data['Vc3'][:,0,0],'+',markersize=2)
plt.plot(U.x,U.data['Vc3'][:,0,0])
plt.title('X field')

plt.figure(5)
#plt.plot(x,u)
plt.plot(V.x,V.data['Vc4'][:,0,0],'+',markersize=2)
plt.plot(U.x,U.data['Vc4'][:,0,0])
plt.title('Y field')


plt.figure(6)
#plt.plot(x,p)
plt.plot(V.x,V.data['Vc5'][:,0,0],'+',markersize=2)
plt.plot(U.x,U.data['Vc5'][:,0,0])
plt.title('Pressure')

plt.ioff()
plt.show()
