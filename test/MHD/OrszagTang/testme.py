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
  inifiles=["idefix.ini",
            "idefix-hll.ini","idefix-hlld-arithmetic.ini",
            "idefix-hlld-hll.ini",
            "idefix-hlld-hlld.ini",
            "idefix-hlld-uct0.ini",
            "idefix-hlld.ini","idefix-tvdlf.ini"]

  for ini in inifiles:
    mytol=tolerance

    test.run(inputFile=ini)
    if test.init and not test.mpi:
      test.makeReference(filename="dump.0001.dmp")

    if(test.single):
      mytol=1e-5

    test.nonRegressionTest(filename="dump.0001.dmp",tolerance=mytol)


test=tst.idfxTest()
if not test.dec:
  test.dec=['2','2']

if not test.all:
  if(test.check):
    test.checkOnly(filename="dump.0001.dmp",tolerance=tolerance)
  else:
    testMe(test)
else:
  for rec in [2,3,4]:
    test.noplot = True
    test.vectPot = False
    test.single=False
    test.reconstruction=rec
    test.mpi=False
    testMe(test)
    test.mpi=True
    testMe(test)

  # single precision validation
  test.reconstruction=2
  test.mpi=False
  test.single=True
  testMe(test)

  # Vector potential validation
  test.single=False
  test.mpi=False
  test.vectPot=True
  testMe(test)
