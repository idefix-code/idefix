#!/usr/bin/env python3
######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/utils/lookupTable/testme.py
#
# Last modified : 08/2026
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG / CNRS - 2023)
# - Clément Robert <cr52@protonmail.com> (2026)
# - Sébastien Valat <sebastien.valat@univ-grenoble-alpes.fr> (IPAG / CNRS - 2026)
# - and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

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
