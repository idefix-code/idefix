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
  inifiles=["idefix-fast.ini","idefix-slow.ini","idefix-alfven.ini","idefix-entropy.ini"]
  for ini in inifiles:
    mytol=tolerance

    test.run(inputFile=ini)
    if test.init and not test.mpi:
      test.makeReference(filename="dump.0001.dmp")
    test.standardTest()
    test.nonRegressionTest(filename="dump.0001.dmp",tolerance=mytol)


test=tst.idfxTest()
if not test.dec:
  test.dec=['2','2','2']

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

  test.mpi=True
  testMe(test)

  test.mpi=False
  test.reconstruction=3
  testMe(test)

  test.reconstruction=4
  testMe(test)
