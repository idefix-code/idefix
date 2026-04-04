#!/usr/bin/env python3
"""

@author: glesur
"""
import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))

import pytools.idfx_test as tst

# Whether we should reset our reference run (only do that on purpose!)

tolerance=1e-13

def testMe(test):
  test.configure()
  test.compile()
  tol=tolerance
  if test.single:
    tol=1e-6

  # default with idefix.ini
  test.run()
  if test.init:
      if not test.mpi:
          test.makeReference(filename="dump.0001.dmp")
  test.nonRegressionTest(filename="dump.0001.dmp",tolerance=tol)


test=tst.idfxTest(__file__)

# if no decomposition is specified, use that one
if not test.dec:
  test.dec=["2","2","2"]

if not test.all:
  if(test.check):
      test.checkOnly(filename="dump.0001.dmp",tolerance=tolerance)
  else:
    testMe(test)
else:
  test.vectPot=False
  test.reconstruction=2
  test.mpi=False
  testMe(test)
  # test in MPI mode
  test.mpi=True
  testMe(test)


  # test with vector potential
  test.mpi=False
  test.vectPot=True
  test.reconstruction=2
  testMe(test)

  test.mpi=True
  testMe(test)

  # test with other precision
  test.single=True
  test.vectPot=False
  test.reconstruction=2
  test.mpi=False
  testMe(test)
  # test in MPI mode
  test.mpi=True
  testMe(test)
