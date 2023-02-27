#!/usr/bin/env python3

"""

@author: glesur
"""
import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))

import pytools.idfx_test as tst

name="dump.0001.dmp"
tolerance=1e-15
def testMe(test):
  test.configure()
  test.compile()
  inifiles=["idefix.ini"]

  # loop on all the ini files for this test
  for ini in inifiles:
    test.run(inputFile=ini)
    test.standardTest();
    if test.init:
      test.makeReference(filename=name)
    test.nonRegressionTest(filename=name,tolerance=tolerance)


test=tst.idfxTest()

if not test.all:
  if(test.check):
    test.checkOnly(filename=name,tolerance=tolerance)
  else:
    testMe(test)
else:
  test.noplot = True

  test.vectPot=False
  test.single=False
  test.reconstruction=2
  test.mpi=False
  testMe(test)
