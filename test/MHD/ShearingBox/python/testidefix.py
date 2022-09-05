#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 15:42:19 2021

@author: lesurg
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import argparse

#compute theoretical solution (from Balbus & Hawley 2006, using notations from Lesur 2021)
def rhs(t, y, Omega, q, B0y, B0z, k0x, k0y, k0z):
    kx = k0x + q*Omega * t * k0y
    ky = k0y
    kz = k0z

    k2 = kx*kx+ky*ky+kz*kz
    vx = y[0]
    vy = y[1]
    vz = y[2]
    bx = y[3]
    by = y[4]
    bz = y[5]

    gxx = kx*kx/k2
    gxy = kx*ky/k2
    gyy = ky*ky/k2
    gyz = ky*kz/k2
    gxz = kx*kz/k2
    #gzz = kz*kz/k2

    kdotB = ky*B0y+kz*B0z

    #we assume b=1j*db hence db=-1j*b

    dvx = kdotB*bx + 2*Omega*vy*(1-gxx) + 2*(1-q)*Omega*gxy*vx
    dvy = kdotB*by - q*Omega*vx*gyy - (2-q)*Omega*(1-gyy)*vx - 2*Omega*vy*gxy
    dvz=  kdotB*bz +2*(1-q)*Omega*vx*gyz-2*Omega*vy*gxz

    dbx = - kdotB*vx
    dby = - kdotB*vy - q*Omega*bx
    dbz = - kdotB*vz

    return np.asarray([dvx,dvy,dvz,dbx,dby,dbz])

parser = argparse.ArgumentParser()
parser.add_argument("-noplot",
                    default=False,
                    help="disable plotting",
                    action="store_true")


args, unknown=parser.parse_known_args()


#initial condition: vr=1, rest is 0, mode initial is nx=0, ny=1, nz=5)
y=solve_ivp(rhs, [0,30], [1, 0, 0, 0, 0, 0], args=(1, 1.5, 0.02, 0.05, 0, 2.0*np.pi, 4*np.pi), dense_output=True)


# read timevol file
rep = "../"
filename = "timevol.dat"

fid=open(rep+filename,"r")

gamma=5/3
# read the first line to get data names
varnames=fid.readline().split()
fid.close()

# load the bulk of the file
data=np.loadtxt(rep+filename,skiprows=1)

# store this in our data structure
V={}

i=0
for name in varnames:
    V[name]=data[:,i]
    i=i+1

# velocity normalisation
v0=V['vx'][0]
# compute L2 error norm
error=np.sqrt((V['vx']/v0-y.sol(V["t"])[0,:])**2 + (V['vy']/v0-y.sol(V["t"])[1,:])**2+ + (V['vz']/v0-y.sol(V["t"])[2,:])**2)



if(not args.noplot):
  # Comparison with exact solution
  plt.rc('text', usetex=True)
  plt.rc('font', family='serif')
  plt.rc('font', size=16)
  plt.close('all')
  plt.figure(1)
  plt.plot(V["t"],V['vx']/v0,'r-',label=r'$u_{R}$')
  plt.plot(V["t"],y.sol(V["t"])[0,:],'r--')
  plt.plot(V["t"],V['vy']/v0,'b-',label=r'$u_{\varphi}$')
  plt.plot(V["t"],y.sol(V["t"])[1,:],'b--')
  plt.plot(V["t"],V['vz']/v0,'g-',label=r'$u_{z}$')
  plt.plot(V["t"],y.sol(V["t"])[2,:],'g--')
  plt.legend()
  plt.xlabel('t')

  # plot error
  plt.figure()
  plt.semilogy(V["t"],error)
  plt.xlim([0,10])
  plt.ylim([1e-3,1])

  plt.ioff()
  plt.show()

err=np.mean(error)
print("Error=",err)

if(err<0.03):
  print("SUCCESS")
  sys.exit(0)
else:
  print("FAILED")
