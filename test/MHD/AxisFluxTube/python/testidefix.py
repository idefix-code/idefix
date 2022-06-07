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

V=readVTKSpherical('../data.0001.vtk')
U=readVTKSpherical('data.0001.ref.vtk')

# Compute BRMS
Brms_ref=np.sqrt(V.data['BX1']**2+V.data['BX2']**2+V.data['BX3']**2)
Brms_sim=np.sqrt(U.data['BX1']**2+U.data['BX2']**2+U.data['BX3']**2)

error=np.mean(np.abs(Brms_ref-Brms_sim)*1e20,axis=(0,1,2))

print("Error=%e"%error)
if error<1e-5:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
