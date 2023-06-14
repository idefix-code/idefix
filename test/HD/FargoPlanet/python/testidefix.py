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
import numpy as np

V=readVTK('../data.0001.vtk', geometry='polar')
U=readVTK('data.0001.ref.vtk', geometry='polar')

# Compute the error on RHO
error=np.mean(np.abs(V.data['RHO']-U.data['rho']),axis=(0,1))

print("Error=%e"%error)
if error<5.0e-2:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
