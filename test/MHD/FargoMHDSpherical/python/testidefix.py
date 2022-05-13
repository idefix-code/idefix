#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:29:41 2020

@author: glesur
"""

import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))
from pytools.vtk_io import readVTKSpherical
import numpy as np

V=readVTKSpherical('../data.0001.vtk')
U=readVTKSpherical('data.0001.ref.vtk')

# Compute the error on BX1
error=np.mean(np.abs(V.data['BX1']-U.data['BX1']))/np.amax(V.data['BX1'])

print("Error=%e"%error)
if error<3.0e-2:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
