#!/usr/bin/env python3

"""

@author: glesur
"""
import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))

import pytools.idfx_test as tst

name="dump.0002.dmp"

tolerance=1e-13

def testMe(test):
  test.configure()
  test.compile()
  inifiles=["idefix.ini"]

  # loop on all the ini files for this test
  for ini in inifiles:
    test.run(inputFile=ini)
    if test.init and not test.mpi:
      test.makeReference(filename=name)
    test.standardTest()
    test.nonRegressionTest(filename=name,tolerance=tolerance)

    # Test the restart option
    dump_mtime = os.path.getmtime("dump.0001.dmp")
    vtk_mtime = os.path.getmtime("data.0005.vtk")

    test.run(inputFile=ini, restart=1)
    test.standardTest()
    test.nonRegressionTest(filename=name,tolerance=tolerance)

    assert dump_mtime == os.path.getmtime("dump.0001.dmp"), "Dump was overwritten on restart"
    assert vtk_mtime == os.path.getmtime("data.0005.vtk"), "VTK file was overwritten on restart"

test=tst.idfxTest()
if not test.dec:
  test.dec=['2','2']

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

  test.mpi=True
  testMe(test)
