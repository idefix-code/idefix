#!/usr/bin/env python3

"""

@author: glesur
"""

import os
import sys

sys.path.append(os.getenv("IDEFIX_DIR"))

import pytools.idfx_test as tst

name = "dump.0001.dmp"

tolerance = 0.0


def testMe(test):
    test.configure()
    test.compile()
    inifiles = [
        "idefix-hll-1.ini",
        "idefix-hll-2.ini",
        "idefix-hllc-1.ini",
        "idefix-hllc-2.ini",
    ]

    # loop on all the ini files for this test
    for ini in inifiles:
        test.run(inputFile=ini)
        if test.init and not test.mpi:
            test.makeReference(filename=name)
        test.standardTest()
        test.nonRegressionTest(filename=name, tolerance=tolerance)


test = tst.idfxTest(__file__)
if not test.dec:
    test.dec = ["2"]

if not test.all:
    if test.check:
        test.checkOnly(filename=name, tolerance=tolerance)
    else:
        testMe(test)
else:
    for rec in range(2, 4):
        test.vectPot = False
        test.single = False
        test.reconstruction = rec
        test.mpi = False
        testMe(test)

    test.mpi = True
    testMe(test)
