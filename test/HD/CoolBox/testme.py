#!/usr/bin/env python3

"""

@author: dutta-alankar
"""
import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))

import pytools.idfx_test as tst

def testMe(test):
  test.configure()
  test.compile()
  test.run(inputFile="idefix.ini")
  test.standardTest()


test=tst.idfxTest(__file__)

if not test.all:
  if(test.check):
    test.standardTest()
  else:
    testMe(test)
else:
  test.noplot = True
  test.single=False
  test.reconstruction=2

  # single core
  test.mpi=False
  testMe(test)

  # MPI
  test.mpi=True
  test.dec=["2","1","2"]
  testMe(test)
