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





# load the dat file produced by the setup
raw=np.loadtxt('../timevol.dat',skiprows=1)
t=raw[:,0]
Ekg=raw[:,1]
Ekd=raw[:,2]
Eth=raw[:,3]

etot=Ekg+Ekd+Eth

if not(args.noplot):
  plt.figure()
  plt.plot(t,Ekg,label="Ek_g")
  plt.plot(t,Ekd,label="Ek_d")
  plt.plot(t,Eth,label="Eth")
  plt.plot(t,etot,'--',label="Etot")
  plt.legend()
  plt.xlabel("t")
  plt.ylabel("Flow energy")
  plt.show()

# Compute relative evolution of total energy
error=abs((etot[-1]-etot[0])/etot[0])

print("error=%f"%error)
if(error<1e-4):
  print("Success!")
else:
  print("Failure!")
  sys.exit(1)
