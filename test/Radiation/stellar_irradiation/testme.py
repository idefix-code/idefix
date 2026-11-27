#!/usr/bin/env python3
######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/Radiation/stellar_irradiation/testme.py
#
# Last modified : 11/2026
#
# Copyright(C) by :
# - Nicolas Scepi <nicolas.scepi@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
# - Sebastien Valat <sebastien.valat@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
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
