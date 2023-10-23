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
  inifiles=["idefix.ini","idefix-rkl.ini","idefix-sl.ini"]
#  inifiles=["idefix-rkl.ini","idefix-sl.ini"]

  # loop on all the ini files for this test
  name="dump.0001.dmp"
  for ini in inifiles:
    test.run(inputFile=ini)
    test.standardTest()
    if test.init:
      test.makeReference(filename=name)
    test.nonRegressionTest(filename=name)


test=tst.idfxTest()

if not test.all:
  if(test.check):
    test.checkOnly(filename="dump.0001.dmp")
  else:
    testMe(test)
else:
  test.noplot = True
  test.vectPot=False
  test.single=False
  test.reconstruction=2
  test.mpi=False
  testMe(test)
