#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:29:41 2020

@author: glesur
"""

import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))
from idefix_pytools.vtk_io import readVTKCart
import numpy as np
import matplotlib.pyplot as plt

V=readVTKCart('../data.0001.vtk')
U=readVTKCart('data.0001.ref.vtk')

# Compute the error on PRS
error=np.mean(np.abs(V.data['PRS']-U.data['prs'])/U.data['prs'],axis=(0,1))

print("Error=%e"%error)
if error<2.5e-2:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
