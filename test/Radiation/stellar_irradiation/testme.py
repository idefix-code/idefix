#!/usr/bin/env python3

"""

@author: glesur
"""

import os
import sys

sys.path.append(os.getenv("IDEFIX_DIR"))

import pytools.idfx_test as tst

name = "dump.0001.dmp"


def testMe(test, mpi):
    inifiles = ["idefix_usertable_kappa.ini"]
    if mpi == True:
        tolerance = 1.0e-10
    else:
        tolerance = 1.0e-13

    # loop on all the ini files for this test
    for ini in inifiles:
        test.configure()
        test.compile()
        test.run(inputFile=ini)
        if test.init and not test.mpi:
            test.makeReference(filename=name)
        test.standardTest()
        test.nonRegressionTest(filename=name, tolerance=tolerance)


test = tst.idfxTest(__file__)

if not test.dec:
    test.dec = ["2", "2"]

if not test.all:
    if test.check:
        test._readLog()
        if test.mpi:
            tolerance = 1.0e-10
        else:
            tolerance = 1.0e-13
        test.checkOnly(filename=name, tolerance=tolerance)
    else:
        testMe(test, test.mpi)
else:
    for rec in [2]:
        test.vectPot = False
        test.single = False
        test.reconstruction = rec
        test.mpi = False
        testMe(test, test.mpi)
