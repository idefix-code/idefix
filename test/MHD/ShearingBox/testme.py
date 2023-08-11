#!/usr/bin/env python3

"""

@author: glesur
"""
import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))

import pytools.idfx_test as tst
tolerance=1e-14

def testMe(test):
  test.configure()
  test.compile()
  inifiles=["idefix.ini","idefix-fargo.ini"]
  for ini in inifiles:
    mytol=tolerance

    test.run(inputFile=ini)
    if test.init and not test.mpi:
      test.makeReference(filename="dump.0001.dmp")
    test.standardTest()
    # When using RKL, except larger error due to B reconstruction and RKL # of substeps
    test.nonRegressionTest(filename="dump.0001.dmp",tolerance=mytol)


test=tst.idfxTest()
if not test.dec:
  test.dec=['2','1','2']

if not test.all:
  if(test.check):
    test.checkOnly(filename="dump.0001.dmp",tolerance=tolerance)
  else:
    testMe(test)
else:
  test.noplot = True
  test.single=False
  test.reconstruction=2
  test.mpi=False
  testMe(test)
