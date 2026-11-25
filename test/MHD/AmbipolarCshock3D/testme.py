#!/usr/bin/env python3
######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/MHD/AmbipolarCshock3D/testme.py
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

tolerance = 3e-14


def testMe(test):
    test.configure()
    test.compile()
    inifiles = ["idefix.ini", "idefix-rkl.ini"]
    for ini in inifiles:
        mytol = tolerance

        test.run(inputFile=ini)
        if test.init and not test.mpi:
            test.makeReference(filename="dump.0001.dmp")
        test.standardTest()
        # When using RKL, except larger error due to B reconstruction and RKL # of substeps
        if test.mpi and ini == "idefix-rkl.ini":
            mytol = 2e-10
        test.nonRegressionTest(filename="dump.0001.dmp", tolerance=mytol)


test = tst.idfxTest(__file__)
if not test.dec:
    test.dec = ["2", "1", "1"]

if not test.all:
    if test.check:
        test.checkOnly(filename="dump.0001.dmp", tolerance=tolerance)
    else:
        testMe(test)
else:
    test.noplot = True
    test.single = False
    test.reconstruction = 2
    test.mpi = False
    testMe(test)

    test.vectPot = True
    testMe(test)

    test.vectPot = False
    test.mpi = True
    testMe(test)
