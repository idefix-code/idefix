#!/usr/bin/env python3

"""

@author: glesur
"""
import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))

import pytools.idfx_test as tst

name="dump.0001.dmp"

def testMe(test):
  test.configure()
  test.compile()
  inifiles=["idefix.ini"]

  # loop on all the ini files for this test
  for ini in inifiles:
    test.run(inputFile=ini)
    if test.init:
      test.makeReference(filename=name)
    test.standardTest()
    test.nonRegressionTest(filename=name,tolerance=1e-14)


test=tst.idfxTest()

if not test.all:
  if(test.check):
    test.checkOnly(filename=name)
  else:
    testMe(test)
else:
  test.noplot = True
  test.reconstruction=2
  testMe(test)
