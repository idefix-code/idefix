#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 15:42:19 2021

@author: lesurg
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as ml


# read timevar file

rep = "../"

fid=open(rep+"/timevol.dat","r")

gamma=5/3
# read the first line to get data names
varnames=fid.readline().split()
fid.close()

# load the bulk of the file
data=np.loadtxt(rep+"/timevol.dat",skiprows=1)

# store this in our data structure
V={}

i=0
for name in varnames:
    V[name]=data[:,i]
    i=i+1

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=16)

 # Enerrgy plots
plt.close('all')
plt.figure(1)
plt.semilogy(V["t"],V['rhovx2'],label=r'$\rho v_x^2$')
plt.semilogy(V["t"],V['rhovy2'],label=r'$\rho v_y^2$')
plt.semilogy(V["t"],V['rhovz2'],label=r'$\rho v_z^2$')
plt.legend()
plt.xlabel('t')

plt.figure(2)
plt.semilogy(V["t"],V['rhovxvy'],label=r'$-\rho v_xv_y$')
plt.semilogy(V["t"],-V['BxBy'],label=r'$B_xB_y$')
plt.legend()
plt.xlabel('t')

W=V['rhovxvy']-V['BxBy']
dt=V['t'][1:]-V['t'][:-1]
theoryP=V['prs'][0]+3/2*(gamma-1)*np.cumsum(dt*W[:-1])
plt.figure(3)
plt.plot(V["t"],V['prs'],label=r'$P$')
plt.plot(V["t"],V['rho'],label=r'$\rho$')
plt.plot(V["t"][:-1],theoryP,label=r'$P_\mathrm{th}$')
plt.legend()
plt.xlabel('t')

plt.figure(4)
plt.plot(V["t"],V['Bx'],label=r'$Bx$')
plt.plot(V["t"],V['By'],label=r'$By$')
plt.plot(V["t"],V['Bz'],label=r'$Bz$')
plt.legend()
plt.xlabel('t')

plt.figure(5)
plt.plot(V["t"],V['Bz']-V['Bz'][0])
plt.xlabel('t')
plt.ylabel(r'$\delta B_z$')

plt.show()
