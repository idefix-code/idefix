#!/usr/bin/env python3

"""

@author: glesur
"""

import os
import sys

sys.path.append(os.getenv("IDEFIX_DIR"))

import pytools.idfx_test as tst

name = "dump.0001.dmp"

tolerance = 1e-10


def testMe(test):
    test.configure()
    test.compile()
    inifiles = [
        "idefix-hll-c001_sub.ini",
        "idefix-hll-c00001_sub.ini",
        "idefix-hllc-c001_sub.ini",
        "idefix-hllc-c00001_sub.ini",
        "idefix-hll-c1_super.ini",
        "idefix-hll-c001_super.ini",
        "idefix-hll-c00001_super.ini",
        "idefix-hllc-c1_super.ini",
        "idefix-hllc-c001_super.ini",
        "idefix-hllc-c00001_super.ini",
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
    for rec in [2]:
        test.vectPot = False
        test.single = False
        test.reconstruction = rec
        test.mpi = False
        testMe(test)
