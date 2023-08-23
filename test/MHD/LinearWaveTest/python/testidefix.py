#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import os
import sys
import math
sys.path.append(os.getenv("IDEFIX_DIR"))
from pytools.dump_io import readDump
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-noplot",
                    default=False,
                    help="disable plotting",
                    action="store_true")


args, unknown=parser.parse_known_args()


def getError(rep):
    try:
        V=readDump(rep+'/dump.0000.dmp')
        U=readDump(rep+'/dump.0001.dmp')
    except:
        return math.nan

    keylist=[]
    keylist=['Vc-RHO','Vc-VX1','Vc-VX2','Vc-VX3','Vs-BX1s','Vs-BX2s','Vs-BX3s','Vc-PRS']

    err = 0.0
    for key in keylist:
        Q1=V.data[key]
        Q2=U.data[key]

        Q1=Q1-np.mean(Q1)
        Q2=Q2-np.mean(Q2)
        print(key+" error=%e"%((np.mean(np.abs(Q1-Q2)))*1e6))
        err=err + (np.mean(np.abs(Q1-Q2)))**2
        if(not args.noplot):
          plt.figure()
          plt.contourf(V.data[key][:,:,16]-U.data[key][:,:,16])
          plt.title(key)


    err = err/len(keylist)
    err = np.sqrt(err)

    return err

error = getError("..")
print("Error=%e"%(error*1e6))

if(not args.noplot):
  plt.show()

if error*1e6<3e-2:
  print("SUCCESS")
  sys.exit(0)
else:
  print("Failed")
  sys.exit(1)
