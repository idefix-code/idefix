#!/usr/bin/env python3

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
  inifiles=["idefix.ini"]

  # loop on all the ini files for this test
  for ini in inifiles:
    test.run(inputFile=ini)
    test.standardTest()
    test.nonRegressionTest(filename="dump.0001.dmp")


test=tst.idfxTest()

if not test.all:
  if(test.check):
    test.checkOnly(filename="dump.0001.dmp")
## The current test is not physics, but merely a non-regression check
## The physical test needs to have the divergence of the velocity field artifiacially set to zero
#  else:
#    testMe(test)
else:
  test.noplot = True
  test.vectPot=False
  test.single=False
  test.reconstruction=2
  test.mpi=False
  testMe(test)
