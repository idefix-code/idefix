import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))

import pytools.idfx_test as tst

test=tst.idfxTest()

rep1="./"
rep2="/Users/lesurg/Documents/src/idefix/test/MHD/AmbipolarCshock3D/"

for n in range(0,1001):
  filename="dump.%04d.dmp"%n
  print("Comparing "+filename)
  test.compareDump(rep1+filename,rep2+filename)
  
