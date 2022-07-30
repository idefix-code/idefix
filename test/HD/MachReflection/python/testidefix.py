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


V=readVTK('../data.0001.vtk', geometry='cartesian')
U=readVTK('data.0001.ref.vtk', geometry='cartesian')

# Compute the error on PRS
error=np.mean(np.abs(np.log(V.data['RHO'])-np.log(U.data['RHO'])),axis=(0,1))

print("Error=%e"%error)
if error<6e-3:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
