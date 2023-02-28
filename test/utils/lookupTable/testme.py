#!/usr/bin/env python3

"""

@author: glesur
"""
import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))
import numpy as np
#from scipy.interpolate import RegularGridInterpolator
import pytools.idfx_test as tst

def MakeNumpyFile():
  x=np.arange(1,10,1.0)
  y=np.arange(5,10,1.0)
  z=np.arange(2,5,1.0)

  xp, yp, zp = np.meshgrid(x,y,z,indexing='ij')

  data=xp+2*yp-zp

  np.save("x.npy",x)
  np.save("y.npy",y)
  np.save("z.npy",z)
  np.save("data.npy",data)
  # show the expected result
  #f=RegularGridInterpolator((x, y, z), data)
  #print(f([2.7,7.4,3.9]))



test=tst.idfxTest()
MakeNumpyFile()

test.configure()
test.compile()
# this test succeeds if it runs successfully
test.run()
