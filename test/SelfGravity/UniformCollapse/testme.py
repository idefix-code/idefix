#!/usr/bin/env python3
######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/SelfGravity/UniformCollapse/testme.py
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

tolerance = 1e-13


def testMe(test):
    test.configure()
    test.compile()
    inifiles = ["idefix.ini"]

    # loop on all the ini files for this test
    for ini in inifiles:
        test.run(inputFile=ini)
        test.standardTest()


test = tst.idfxTest(__file__)

if not test.all:
    if test.check:
        test.standardTest()
    else:
        testMe(test)
else:
    test.noplot = True
    test.vectPot = False
    test.single = False
    test.reconstruction = 2
    test.mpi = False
    testMe(test)
