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

  test.run(inputFile="idefix.ini")
  if test.init and not test.mpi:
    test.makeReference(filename="dump.0001.dmp")

  test.nonRegressionTest(filename="dump.0001.dmp",tolerance=tolerance)


test=tst.idfxTest()
if not test.dec:
  test.dec=['2','2']

if not test.all:
  if(test.check):
    test.checkOnly(filename="dump.0001.dmp",tolerance=tolerance)
  else:
    testMe(test)
else:
  test.noplot = True
  test.vectPot = False
  test.single=False
  test.reconstruction=2
  test.mpi=False
  testMe(test)
  test.mpi=True
  testMe(test)
