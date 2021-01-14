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
from idefix_testing.framework import readVTKPolar
import numpy as np
import matplotlib.pyplot as plt

V=readVTKPolar('../data.0001.vtk')
U=readVTKPolar('data.0001.ref.vtk')

# Compute the error on PRS
error=np.mean(np.abs(V.data['VX2']-U.data['vx2']),axis=(0,1))

print("Error=%e"%error)
if error<1.5e-3:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
