#!/usr/bin/env python3

"""

@author: glesur
"""
import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))

import pytools.idfx_test as tst

name="dump.0001.dmp"

test=tst.idfxTest()

if not test.dec:
  test.dec=['2','2','2']

if test.check:
  print("nothing to do for the check, just check that the xdmf file exists")
else:
  test.vectPot=False
  test.single=False
  test.reconstruction=2
  test.mpi=False
  # Only check that the test runs.
  test.configure(definitionFile="definitions.hpp")
  test.compile()
  test.run(inputFile="idefix.ini")
