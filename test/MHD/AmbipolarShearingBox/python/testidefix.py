#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 15:31:57 2020

@author: lesurg
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))
from pytools.vtk_io import readVTK

rep='../'
nend=1000

n=0
dt=0.1


Bx=np.zeros(nend)

t=dt*np.arange(0,nend)
for n in range(nend):
    V=readVTK(rep+'/data.'+'%0*d'%(4,n)+'.vtk', geometry='cartesian')
    Bx[n]=np.sqrt(np.mean(np.mean(V.data['BX1']**2,axis=2),axis=0))

plt.figure()
plt.semilogy(t,Bx)
plt.semilogy(t,Bx[-1]*np.exp(0.171*(t-t[-1])))

gr=np.log(Bx[-1]/Bx[-700])/(t[-1]-t[-700])
print('growth rate=%g'%gr)
if np.abs(gr-0.171)/0.171<0.05:
    print("SUCCESS!")
else:
    print("FAILED!")
