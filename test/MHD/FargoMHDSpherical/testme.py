#!/usr/bin/env python3

"""

@author: glesur
"""
import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))

import pytools.idfx_test as tst

name="dump.0001.dmp"

tolerance=1e-14

def testMe(test):
  test.configure()
  test.compile()
  inifiles=["idefix.ini"]

  # loop on all the ini files for this test
  for ini in inifiles:
    test.run(inputFile=ini)
    if test.init and not test.mpi:
      test.makeReference(filename=name)
    test.nonRegressionTest(filename=name,tolerance=tolerance)


test=tst.idfxTest()
if not test.dec:
  test.dec=['2','2','2']

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

  #test vector potential formulation
  test.vectPot=True
  testMe(test)

  # test with MPI
  test.vectPot=False
  test.mpi=True
  testMe(test)
