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
import argparse
import inifix

parser = argparse.ArgumentParser()
parser.add_argument("-i",
                    default=False,
                    help="idefix input file")

args, unknown=parser.parse_known_args()

referenceFile = 'data.0001.ref.vtk'

if(args.i):
  conf = inifix.load(args.i)
  if "coarsening" in conf["Grid"]:
    referenceFile = 'data.0001.ref-coarsening.vtk'

V=readVTK('../data.0001.vtk')
U=readVTK(referenceFile)

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
