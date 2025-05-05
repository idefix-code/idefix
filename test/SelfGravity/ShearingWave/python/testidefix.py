#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:42:19 2025

@author: vdbma
"""
import os
import sys
import numpy as np
import inifix
sys.path.append(os.getenv("IDEFIX_DIR"))
from pytools.vtk_io import readVTK
import argparse
import scipy.integrate as si

# check that the potential minima is sheared
# at the correct rate in the boundaries


parser = argparse.ArgumentParser()
parser.add_argument("-noplot",
                    default=False,
                    help="disable plotting",
                    action="store_true")


args, unknown=parser.parse_known_args()

conf = inifix.load('../idefix.ini')

# unit dependant stuff
Omega = conf["Hydro"]["rotation"]
G = conf["Gravity"]["gravCst"]
# user parameters
q = -conf["Hydro"]["shearingBox"]
isSG = "potential" in conf["Gravity"]


kappa = 2*Omega*Omega*(1-q)
sigma0 = 0.0005
Sigma0 = 1/40
kx0 = -4*np.pi
ky = 2*np.pi
Q = 1


xi0 = (2-q)*Omega/Sigma0
xi1 = -sigma0*xi0
tmax = conf["TimeIntegrator"]["tstop"]

cs = np.pi*G*Sigma0*Q/kappa


def f2s(t,y):
    kx = kx0 + q*Omega*ky*t
    k = np.sqrt(kx*kx+ky*ky)

    pot =  -4*np.pi *G *Sigma0 *isSG

    b = -2*q*Omega*kx*ky/k/k
    c = 2*(2-q)*Omega*Omega+k**2*cs**2 + pot - 2*q*(2-q)*Omega*Omega*ky*ky/k/k
    d = 2*Omega*(1-q*ky*ky/k/k)*Sigma0*xi1

    sp,s = y
    f1 = -b*sp-c*s-d
    f2 = sp
    return np.array([f1,f2])


y0 = [0,sigma0]


amplitude = []
amp_exact = []
time = []

sol = si.solve_ivp(f2s,y0=y0,t_span = [0,tmax],t_eval=np.linspace(0,tmax,2048))


for i in range(200):
    V = readVTK(f"../data.{i:04}.vtk",geometry="cartesian")
    t=V.t[0]
    xx,yy = np.meshgrid(V.x,V.y,indexing="ij")
    time.append(t)

    mode = np.exp(1j*((kx0+q*Omega*t*ky)*xx + ky*yy))
    amplitude.append(np.sum(V.data["RHO"][:,:,0]*mode.conjugate())/V.nx/V.ny*2) # 2 parce que int cos**2 = 1/2
    i_t = np.argmin(np.abs(t-sol.t))
    amp_exact.append(sol.y[1][i_t])


time = np.array(time)
amplitude=np.array(amplitude).real/Sigma0
amp_exact=np.array(amp_exact)



err=np.sum((amp_exact-amplitude)**2)/np.sum(amp_exact**2)
print("Error=",err)


if not args.noplot:
  import matplotlib.pyplot as plt
  fig, ax = plt.subplots()
  ax.plot(sol.t, sol.y[1],color="orange",ls="-",label="linear")
  ax.plot(time,amplitude,".",label="vtk")
  ax.legend()
  plt.show()

if(err<.1):
  print("SUCCESS")
  sys.exit(0)
else:
  print("FAILED")
  sys.exit(1)
