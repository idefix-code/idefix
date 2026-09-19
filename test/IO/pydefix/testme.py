#!/usr/bin/env python3
######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/IO/pydefix/testme.py
#
# Last modified : 07/2026
#
# Copyright(C) by :
# - Clément Robert <clement.robert@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
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

tolerance = 1e-12


def testMe(test):
    test.configure()
    test.compile()

    test.run(inputFile="idefix.ini")
    if test.init and not test.mpi:
        test.makeReference(filename="dump.0001.dmp")

    test.nonRegressionTest(filename="dump.0001.dmp", tolerance=tolerance)


test = tst.idfxTest(__file__)
if not test.dec:
    test.dec = ["2", "2"]

if not test.all:
    if test.check:
        test.checkOnly(filename="dump.0001.dmp", tolerance=tolerance)
    else:
        testMe(test)
else:
    test.noplot = True
    test.vectPot = False
    test.single = False
    test.reconstruction = 2
    test.mpi = False
    testMe(test)
    test.mpi = True
    testMe(test)
