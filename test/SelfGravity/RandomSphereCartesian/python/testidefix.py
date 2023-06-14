#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:29:41 2020

@author: glesur
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
import argparse

sys.path.append(os.getenv("IDEFIX_DIR"))

from pytools.vtk_io import readVTK

# Compute theoretical psi from spectral solution to Laplace equation

def getPsiTheoretical(x,y,z,rho):

  n=x.shape[0]

  # Compute fourier transform and k vector
  k1D=2.0*np.pi*fft.fftfreq(n,1/n)
  [kx,ky,kz]=np.meshgrid(k1D,k1D,k1D,indexing='ij')
  k2=(kx**2+ky**2+kz**2)
  # avoid divide by 0
  k2[0,0,0]=1.0
  ik2=1.0/k2
  # kill zero frequency
  ik2[0,0,0]=0

  rhof=fft.fft2(rho,axes=(0,1,2))
  psif=-rhof*ik2
  psi=4.0*np.pi*np.real(fft.ifft2(psif,axes=(0,1,2)))
  return(psi)

# Add noplot argument to parser
parser = argparse.ArgumentParser()
parser.add_argument("-noplot",
                    default=False,
                    help="disable plotting",
                    action="store_true")

args=parser.parse_args()

# Read vtk
V=readVTK('../data.0000.vtk')



# Cartesian coordinates
x = V.x
y = V.y
z = V.z

# Density (for norm)
rho = V.data['RHO']

# Constrained potential
phi_f = V.data['phiP']

phi_th=getPsiTheoretical(x,y,z,rho)



# Recovering gradients to get rid of the integration constant
gx, gy, gz = np.zeros(phi_f.shape), np.zeros(phi_f.shape), np.zeros(phi_f.shape)
gx_f=-np.gradient(phi_f,axis=0)
gx_th=-np.gradient(phi_th,axis=0)
gy_f=-np.gradient(phi_f,axis=1)
gy_th=-np.gradient(phi_th,axis=1)
gz_f=-np.gradient(phi_f,axis=2)
gz_th=-np.gradient(phi_th,axis=2)

norm=np.sqrt(np.mean((gx_f-gx_th)**2+(gy_f-gy_th)**2+(gz_f-gz_th)**2)/np.mean(gx_th**2+gy_th**2+gz_th**2))

# Display a check plot if necessary
if(not args.noplot):
    n=x.size//2
    plt.plot(x,gx_f[:,n,n], 'c-', label='cpt')
    plt.plot(x, gx_th[:,n,n], 'k--', label='th')
    plt.xlabel('r')
    plt.ylabel('g')
    plt.legend()
    plt.show()

# Print the result of the test
print("Error=%e"%norm)
if norm<1e-2:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
