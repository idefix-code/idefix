#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 09:10:06 2020

@author: lesurg
"""
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))
from pytools.vtk_io import readVTK
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-noplot",
                    default=False,
                    help="disable plotting",
                    action="store_true")


args, unknown=parser.parse_known_args()

# Parameters from Maclow+ 1995
theta=np.pi/4
M=50
b0=np.sin(theta)
A=10
L=1

xcut=40.0

def f(x,D):
    b=np.sqrt(b0**2+2*A**2*(D-1)*(1/D-1/M**2))
    diff=(b-b0)/A**2*np.cos(theta)**2+np.sin(theta)
    diff=(b-D*diff)/(b**2+np.cos(theta)**2)
    diff=(b/A)*diff/(L*(1/D**2-1/M**2))

    return(diff)

# load solution
V=readVTK('../data.0001.vtk', geometry='cartesian')


if(V.y.size>V.x.size):
    x=V.y
    DSim=V.data['RHO'][0,:,0]
    bxSim=V.data['BX2'][0,:,0]
    bySim=V.data['BX3'][0,:,0]
elif(V.z.size>V.x.size):
    x=V.z
    DSim=V.data['RHO'][0,0,:]
    bxSim=V.data['BX3'][0,0,:]
    bySim=V.data['BX1'][0,0,:]
else:
    x=V.x
    DSim=V.data['RHO'][:,0,0]
    bxSim=V.data['BX1'][:,0,0]
    bySim=V.data['BX2'][:,0,0]

bSim=bySim/np.sqrt(bySim[0]**2+bxSim[0]**2)

# Index where the two solutions are assumed to match
iref=np.argwhere(DSim>5.0)[0][0]

#spatial index where we cut the solution
iend=np.argwhere(x>xcut)[0][0]

# compute analytical solution
r=ode(f).set_integrator('vode',rtol=1e-6)
r.set_initial_value(DSim[iref],x[iref])

Dth=np.zeros(x.shape)
for i in range(x.size):
    r.set_initial_value(DSim[iref],x[iref])
    Dth[i]=r.integrate(x[i]).item()
    #print("Dth[%d]=%g"%(i,Dth[i]))

bTh=np.sqrt(b0**2+2*A**2*(Dth-1)*(1/Dth-1/M**2))

errb=(bSim[:iend]-bTh[:iend])/np.amax(bTh[:iend])
errD=(DSim[:iend]-Dth[:iend])/np.amax(Dth[:iend])

if(not args.noplot):
    plt.close('all')
    plt.figure()
    plt.subplot(211)
    plt.plot(x[:iend]/L,bSim[:iend],'--',label='Sim')
    plt.plot(x[:iend]/L,bTh[:iend],label='Theoretical')
    plt.xlabel('x/L')
    plt.ylabel('b')
    plt.subplot(212)
    plt.plot(x[:iend]/L,errb,'--',label='Error')
    plt.xlabel('x/L')
    plt.ylabel('b')
    plt.legend()

    plt.figure()
    plt.subplot(211)
    plt.plot(x[:iend]/L,DSim[:iend],label='Sim')
    plt.plot(x[:iend]/L,Dth[:iend],'--',label='Theoretical')
    plt.xlabel('x/L')
    plt.ylabel('D')
    plt.subplot(212)
    plt.plot(x[:iend]/L,errD,'--',label='Error')
    plt.xlabel('x/L')
    plt.ylabel('D')
    plt.legend()

    plt.ioff()
    plt.show()

err=np.mean(0.5*np.sqrt(errb**2+errD**2))
print("Error total=%e"%err)
if(err<5.9e-3):
    print("Success!")
    sys.exit(0)
else:
    print("Failed!")
    sys.exit(1)
