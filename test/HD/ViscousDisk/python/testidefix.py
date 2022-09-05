#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:29:41 2020

@author: glesur
"""

import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))
from pytools.vtk_io import readVTK
import numpy as np
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-noplot",
                    default=False,
                    help="disable plotting",
                    action="store_true")


args, unknown=parser.parse_known_args()


V=readVTK('../data.0001.vtk', geometry='spherical')


iref=32
alpha=0.002
epsilon=0.1
# q and p following Fromang+2011
q=-1
p=-1.5


vr=V.data['VX1'][iref,:,0]
vphi=V.data['VX3'][iref,:,0]
r=V.r[iref]
R=V.r[iref]*np.sin(V.theta)
Z=V.r[iref]*np.cos(V.theta)

Zh=Z/(epsilon*R)

# Theoretical profile
ur=-alpha*epsilon**2*R**(q+1/2)*(3*p+2*q+6+(5*q+9)/2*Zh**2)
#uphi=R**(-0.5)*(1+1/2*epsilon**2*(p+q+q/2*Zh**2))
uphi=R**(-0.5)*(R/r-2.5*epsilon*epsilon)**0.5

# Compute error on accretion speed
error=np.mean(np.fabs(ur-vr))/np.max(np.fabs(ur))

if(not args.noplot):
    plt.figure()
    plt.plot(Zh,vr,label='run')
    plt.plot(Zh,ur,label='theoretical')
    plt.xlabel('Z/h')
    plt.ylabel('Vr')
    plt.legend()

    plt.figure()
    plt.plot(Zh,vphi,label='run')
    plt.plot(Zh,uphi,label='theoretical')
    plt.xlabel('Z/h')
    plt.ylabel('Vphi')
    plt.legend()

    # check radial equibrium
    n=V.theta.size//2
    rhomid=V.data['RHO'][:,n,0]
    R=V.r*np.sin(V.theta[n])
    cs=epsilon*R**(-0.5)
    Pmid=rhomid*cs*cs
    Vmid=V.data['VX3'][:,n,0]
    dV2=(Vmid**2-V.r**(-1))/R*rhomid
    dP=(Pmid[2:]-Pmid[:-2])/(V.r[2]-V.r[0])

    plt.figure()
    plt.plot(V.r[1:-1],dP,label="$d_r P$")
    plt.plot(V.r,dV2,label="$(V^2-V_K^2)/R$")
    plt.legend()
    plt.xlabel('R')
    plt.ylabel('Radial forces in the midplane')

    plt.show()
#plt.figure()
#plt.plot(V.r,rhomid,V.r,R**(-1.5))



print("Error=%e"%error)
if error<1.1e-2:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
