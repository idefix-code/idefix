#!/usr/bin/env python3
######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/Planet/Planet3Body/testme.py
#
# Last modified : 07/2026
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG / CNRS - 2023)
# - Clément Robert <cr52@protonmail.com> (2026)
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

import pytools.idfx_test as tst

name = "dump.0001.dmp"

tolerance = 2e-11


def testMe(test):
    test.configure()
    test.compile()
    inifiles = ["idefix.ini"]

    # loop on all the ini files for this test
    for ini in inifiles:
        test.run(inputFile=ini)
        if test.init and not test.mpi:
            test.makeReference(filename=name)
        test.standardTest()
        test.nonRegressionTest(filename=name, tolerance=tolerance)


test = tst.idfxTest(__file__)

if not test.all:
    if test.check:
        test.checkOnly(filename=name, tolerance=tolerance)
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
