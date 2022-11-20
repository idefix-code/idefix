#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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


# values from the setup to compute the decay rate
kappa=0.1
gamma=1.4
rho0=1.0

#diffusion coefficient
D=kappa*(gamma-1)/rho0

k=2.0*np.pi

# decay rate of the temperature perturbation
gamma=-k**2*D

raw=np.loadtxt('../analysis.dat',skiprows=1)
t=raw[:,0]
Tidfx=raw[:,1]

Tth=Tidfx[0]*np.exp(gamma*t)
error=np.mean(Tidfx/Tth-1.0)

if(not args.noplot):

    plt.close('all')
    plt.figure()
    plt.plot(t,Tidfx/Tth-1.0)

    print("Error=%e"%error)
    plt.ioff()
    plt.show()

if(error<2.5e-6):
    print("SUCCESS")
    sys.exit(0)
else:
    print("Failed")
    sys.exit(1)
