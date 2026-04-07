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

  # Check restarts
  test.run("idefix.ini")


test=tst.idfxTest()

# if no decomposition is specified, use that one
if not test.dec:
  test.dec=["2","2","2"]

if not test.all:
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
