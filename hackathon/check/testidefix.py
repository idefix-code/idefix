#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:29:41 2020

@author: glesur
"""

import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))
from pytools.vtk_io import readVTKPolar
import numpy as np
import matplotlib.pyplot as plt

print("Checking results consistency");

V=readVTKPolar('../data.0001.vtk')
U=readVTKPolar('data.ref.vtk')

# Compare fields
errormax=0

for key in U.data.keys():
    error=np.mean(np.abs(V.data[key]-U.data[key]),axis=(0,1,2))
    print("Error on "+key+"=%e"%error)
    if error>errormax:
        errormax=error

if(errormax<1e-13):
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
