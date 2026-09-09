#!/usr/bin/env python3
######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/MHD/sod-iso/testme.py
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


def testMe(test):
    test.configure()
    test.compile()
    inifiles = ["idefix.ini", "idefix-hll.ini", "idefix-hlld.ini", "idefix-tvdlf.ini"]
    if test.reconstruction == 4:
        inifiles = ["idefix-rk3.ini", "idefix-hlld-rk3.ini"]

    # loop on all the ini files for this test
    for ini in inifiles:
        test.run(inputFile=ini)
        if test.init:
            test.makeReference(filename=name)
        test.nonRegressionTest(filename=name)


test = tst.idfxTest(__file__)

if not test.all:
    if test.check:
        test.checkOnly(filename=name)
    else:
        testMe(test)
else:
    test.noplot = True
    for rec in range(2, 5):
        test.vectPot = False
        test.single = False
        test.reconstruction = rec
        test.mpi = False
        testMe(test)

    # test in single precision
    test.reconstruction = 2
    test.single = True
    testMe(test)
