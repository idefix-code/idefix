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
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

parser = argparse.ArgumentParser()
parser.add_argument("-noplot",
                    default=False,
                    help="disable plotting",
                    action="store_true")


args, unknown=parser.parse_known_args()

V=readVTK('../data.0001.vtk', geometry='cartesian')

# Find where the shock starts
istart = np.argwhere(V.data['RHO'][:,0,0]>2.0)[0][0]-1

# Compute analytical solution
M=2 # Mach number
K=np.asarray([1,3,5])/M # Drag coefficients

# (56) from Benitez-Llambay 2018
def get_wg(Mach, wi):
    coeff=[1,np.sum(wi-1)-1/Mach**2-1,1/Mach**2]
    solution=np.roots(coeff)
    return(np.min(solution))

# (55) from Benitez-Llambay 2018
def rhs(t,wi,Mach,K):
    wg = get_wg(Mach,wi)
    return K*(wg-wi)

wi0=np.asarray([1.0,1.0,1.0])
x=V.x[istart:]
arguments=M,K

sol = solve_ivp(rhs, t_span=[x[0],x[-1]], y0=wi0, t_eval=x,args=arguments)

# compute gas velocity
wg=np.zeros(len(x))
for i in range(0,len(wg)):
    wg[i]=get_wg(M,sol.y[:,i])


if(not args.noplot):
    plt.figure(1)
    plt.plot(V.x,V.data['RHO'][:,0,0],'o',markersize=4,color='tab:purple')
    plt.plot(x,1/wg,color='tab:purple')
    plt.plot(V.x,V.data['Dust0_RHO'][:,0,0],'o',markersize=4,color='tab:orange')
    plt.plot(x,1/sol.y[0,:],color='tab:orange')
    plt.plot(V.x,V.data['Dust1_RHO'][:,0,0],'o',markersize=4,color='tab:blue')
    plt.plot(x,1/sol.y[1,:],color='tab:blue')
    plt.plot(V.x,V.data['Dust2_RHO'][:,0,0],'o',markersize=4,color='tab:green')
    plt.plot(x,1/sol.y[2,:],color='tab:green')
    plt.xlim([0,10])
    plt.title('Density')

    plt.figure(2)
    plt.plot(V.x,V.data['VX1'][:,0,0],'o',markersize=4,color='tab:purple')
    plt.plot(x,M*wg,color='tab:purple')
    plt.plot(V.x[::4],V.data['Dust0_VX1'][::4,0,0],'o',markersize=4,color='tab:orange')
    plt.plot(x,M*sol.y[0,:],color='tab:orange')
    plt.plot(V.x[::4],V.data['Dust1_VX1'][::4,0,0],'o',markersize=4,color='tab:blue')
    plt.plot(x,M*sol.y[1,:],color='tab:blue')
    plt.plot(V.x[::4],V.data['Dust2_VX1'][::4,0,0],'o',markersize=4,color='tab:green')
    plt.plot(x,M*sol.y[2,:],color='tab:green')

    plt.xlim([0,10])
    plt.title('Velocity')

    plt.ioff()
    plt.show()

# Compute the error on the dust densities

error=np.mean(np.fabs(V.data['Dust0_RHO'][istart:,0,0]-1/sol.y[0,:]) + np.fabs(V.data['Dust1_RHO'][istart:,0,0]-1/sol.y[1,:]) + np.fabs(V.data['Dust2_RHO'][istart:,0,0]-1/sol.y[2,:])) / 3
print("Error=%e"%error)
if error<3.6e-2:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
