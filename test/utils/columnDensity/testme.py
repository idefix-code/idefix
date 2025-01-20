#!/usr/bin/env python3

"""

@author: glesur
"""
import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))
import pytools.idfx_test as tst

test=tst.idfxTest()

test.configure()
test.compile()
# this test succeeds if it runs successfully
test.run()

test.mpi = True
test.configure()
test.compile()
# this test succeeds if it runs successfully
test.run()
