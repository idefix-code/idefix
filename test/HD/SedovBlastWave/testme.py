#!/usr/bin/env python3
######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/HD/SedovBlastWave/testme.py
#
# Last modified : 07/2026
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2023)
# - Clément Robert <clement.robert@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

"""

@author: glesur
"""

import os
import sys

sys.path.append(os.getenv("IDEFIX_DIR"))

import pytools.idfx_test as tst

name = "dump.0001.dmp"

test = tst.idfxTest(__file__)

if not test.dec:
    test.dec = ["2", "2", "2"]

if test.check:
    test.standardTest()
else:
    test.vectPot = False
    test.single = False
    test.reconstruction = 2
    test.mpi = True
    # Cartesian validation
    test.configure(definitionFile="definitions.hpp")
    test.compile()
    test.run(inputFile="idefix.ini")
    test.standardTest()

    # Spherical validation
    test.configure(definitionFile="definitions-spherical.hpp")
    test.compile()
    test.run(inputFile="idefix-spherical.ini")
    test.standardTest()
