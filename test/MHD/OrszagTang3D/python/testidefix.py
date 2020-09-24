#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:29:41 2020

@author: glesur
"""

import idefixTools as idfx
import numpy as np
import matplotlib.pyplot as plt

V=idfx.readVTKCart('../data.0005.vtk')
U=idfx.readVTKCart('data.0005.ref.vtk')

# Compute the error on PRS
error=np.sqrt(np.mean((V.data['PRS']-U.data['PRS'])**2/V.data['PRS']**2,axis=(0,1,2)))

print("Error=%e"%error)
if error<1e-2:
    print("SUCCESS!")
else:
    print("FAILURE!")
