#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:29:41 2020

@author: glesur
"""

import os
import sys
TESTDIR_PATH = os.path.join(os.getenv("IDEFIX_DIR"), "test")
sys.path.append(TESTDIR_PATH)
from idefix_testing.framework import readVTKCart
import numpy as np
import matplotlib.pyplot as plt

V=readVTKCart('../data.0001.vtk')
U=readVTKCart('data.0001.ref.vtk')

# Compute the error on PRS
error=np.sqrt(np.mean((V.data['PRS']-U.data['PRS'])**2/V.data['PRS']**2,axis=(0,1,2)))

print("Error=%e"%error)
if error<3e-3:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
