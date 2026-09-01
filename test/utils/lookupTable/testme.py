#!/usr/bin/env python3

"""

@author: glesur
"""

import os
import sys

sys.path.append(os.getenv("IDEFIX_DIR"))
# from scipy.interpolate import RegularGridInterpolator
import testmelib

import pytools.idfx_test as tst

test = tst.idfxTest(__file__)
testmelib.MakeNumpyFile()

test.configure()
test.compile()
# this test succeeds if it runs successfully
test.run()
