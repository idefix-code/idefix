#!/usr/bin/env python3

"""

@author: glesur
"""
import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))

import pytools.idfx_test as tst
tolerance=1e-12
def testMe(test):
  test.configure()
  test.compile()
  inifiles=["idefix.ini"]#,"idefix-noSG.ini"]

  for ini in inifiles:
    test.run(inputFile=ini)
    if test.init and not test.mpi:
      test.makeReference(filename="dump.0001.dmp")
    test.standardTest()
    test.nonRegressionTest(filename="dump.0001.dmp",tolerance=tolerance)


test=tst.idfxTest()
test.noplot = False

if not test.dec:
  test.dec=['8','1']

if not test.all:
  if(test.check):
    test.checkOnly(filename="dump.0001.dmp",tolerance=tolerance)
  else:
    testMe(test)
else:
  test.single=False
  test.reconstruction=2
  test.mpi=False
  testMe(test)

  test.mpi=True
  testMe(test)
