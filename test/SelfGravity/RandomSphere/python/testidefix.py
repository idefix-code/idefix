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
import argparse

sys.path.append(os.getenv("IDEFIX_DIR"))

from pytools.vtk_io import readVTK

# Add noplot argument to parser
parser = argparse.ArgumentParser()
parser.add_argument("-noplot",
                    default=False,
                    help="disable plotting",
                    action="store_true")

args=parser.parse_args()

# Read vtk
V=readVTK('../data.0000.vtk')

# Spherical coordinates
nr, nth, nphi = np.squeeze(V.nx), np.squeeze(V.ny), np.squeeze(V.nz)
r = np.squeeze(V.r)
th = np.squeeze(V.theta)
phi = np.squeeze(V.phi)
rg, thg, phig = np.meshgrid(r, th, phi, indexing="ij")

# Cartesian coordinates
x = rg*np.sin(thg)*np.cos(phig)
y = rg*np.sin(thg)*np.sin(phig)
z = rg*np.cos(thg)

# Density (for norm)
rho = np.squeeze(V.data['RHO'])

# Constrained potential
phi_f = np.squeeze(V.data['phiP'])

# Theoretical parameters
r0 = 3.
delta = 1.
rho0 = 3./(4.*np.pi)
G = 1.
M = 4.*np.pi/3.*delta**3*rho0

# Theoretical potential
dr = np.sqrt((x-r0)**2 + y**2 + z**2)
sph = np.where(dr <= delta) # Sphere delimitation
phi_th = -G*M/dr # Out of the sphere
phi_th[sph] = - (G*M)/(2*delta**3)*(3*delta**2 - dr[sph]**2) # In the sphere

# Recovering gradients to get rid of the integration constant
gth, gf = np.zeros(phi_f.shape), np.zeros(phi_f.shape)
for i in range(nphi):
    for j in range(nth):
        gth[:,j,i] = -np.gradient(phi_th[:,j,i], r, axis=0)
        gf[:,j,i] = -np.gradient(phi_f[:,j,i], r, axis=0)

# Recovering the norm (density ponderated like the solver)
num = np.linalg.norm(gf-gth)
denom = np.linalg.norm(rho)
norm = num/denom
print(norm)

# Display a check plot if necessary
if(not args.noplot):
    plt.plot(r, gf[:,nth//2-1,0], 'c-', label='cpt')
    plt.plot(r, gth[:,nth//2-1,0], 'k--', label='th')
    plt.title('Gravitational gradient along the sphere solid angle')
    plt.xlabel('r')
    plt.ylabel('g')
    plt.legend()
    plt.show()

# Print the result of the test
print("Error=%e"%norm)
if norm<1.5e-1:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
