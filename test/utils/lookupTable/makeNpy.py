#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:29:41 2020

@author: glesur
"""


import numpy as np
from scipy.interpolate import RegularGridInterpolator

x=np.arange(1,10,1.0)
y=np.arange(5,10,1.0)
z=np.arange(2,5,1.0)

xp, yp, zp = np.meshgrid(x,y,z,indexing='ij')

data=xp+2*yp-zp

np.save("x.npy",x)
np.save("y.npy",y)
np.save("z.npy",z)
np.save("data.npy",data)
f=RegularGridInterpolator((x, y, z), data)
print(f([2.7,7.4,3.9]))
