#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:29:41 2020

@author: glesur
"""

import idefixTools as idfx
import numpy as np
import sys
import matplotlib.pyplot as plt

V=idfx.readVTKCart('../data.0001.vtk')
U=idfx.readVTKCart('data.0001.ref.vtk')

# Compute the error on PRS
error=np.mean(np.abs(np.log(V.data['RHO'])-np.log(U.data['RHO'])),axis=(0,1))

print("Error=%e"%error)
if error<6e-3:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
