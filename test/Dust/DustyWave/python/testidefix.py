#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 14:45:58 2023

@author: lesurg
"""
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-noplot",
                    default=False,
                    help="disable plotting",
                    action="store_true")


args, unknown=parser.parse_known_args()


#  solve eq. A6 in Riols  & Lesur (2018)
cs=1.0
chi=1.0
kx=2.0*np.pi
taus=1.0
poly=np.zeros(4,dtype=complex)
poly[0]=1
poly[1]=1j/taus*(1+chi)
poly[2]=-kx**2*cs**2
poly[3]=-1j*kx**2*cs**2/taus
sol=np.roots(poly)

# Get the minimal decay rate (this should be the one that pops up)
tau=np.amax(np.imag(sol))




# load the dat file produced by the setup
raw=np.loadtxt('../timevol.dat',skiprows=1)
t=raw[:,0]
vx2=raw[:,1]
rho2=raw[:,2]

etot=vx2+rho2

if not(args.noplot):
  plt.figure()
  plt.semilogy(t,etot,label="idefix")
  plt.semilogy(t,np.exp(2*tau*t)*etot[10],'--',label="theoretical decay rate")
  plt.legend()
  plt.xlabel("t")
  plt.ylabel("Wave energy")
  plt.show()

# Compute decay rate:
tau_measured=t[-1]/(2*np.log(etot[-1]/etot[0]))
# error on the decay rate:
error=(tau_measured-tau)/tau

print("error=%f"%error)
if(error<0.065):
  print("Success!")
else:
  print("Failure!")
  sys.exit(1)
