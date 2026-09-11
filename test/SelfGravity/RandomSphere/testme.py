#!/usr/bin/env python3
######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/SelfGravity/RandomSphere/testme.py
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


def testMe(test):
    test.configure()
    test.compile()
    inifiles = ["idefix.ini", "idefix-cg.ini", "idefix-minres.ini"]

    # loop on all the ini files for this test
    for ini in inifiles:
        test.run(inputFile=ini)
        # if test.init:
        #  test.makeReference(filename=name)
        test.standardTest()
        # since the gravitationnal potential is not included in .dmp files, we can't perform
        # the full non-regression test
        # test.nonRegressionTest(filename=name)


test = tst.idfxTest(__file__)

if not test.dec:
    test.dec = ["2", "2", "1"]

if not test.all:
    testMe(test)
else:
    test.mpi = False
    test.noplot = True
    testMe(test)
    test.mpi = True
    testMe(test)
