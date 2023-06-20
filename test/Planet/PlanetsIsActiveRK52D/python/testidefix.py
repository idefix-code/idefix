#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  21 12:04:41 2021

@author: gwafflard
"""

import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))
import numpy as np

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def datafile(filename):
    return np.loadtxt(filename, dtype="float64").T

def time_dist(filename):
    columns = datafile(filename)
    return(columns[-1], np.sqrt(columns[1]**2+columns[2]**2+columns[3]**2))

plot = False

time0_active, active0 = datafile("../isactive0.dat")
time1_active, active1 = datafile("../isactive1.dat")

i0_active = find_nearest(time0_active, 0.0)
i1_active = find_nearest(time1_active, 15.0)

time0_rk4, dist0_rk4 = time_dist("planet0.rk4.dat")
time0_rk5, dist0_rk5 = time_dist("planet0.rk5.dat")
time0, dist0 = time_dist("../planet0.dat")
time1_rk4, dist1_rk4 = time_dist("planet1.rk4.dat")
time1_rk5, dist1_rk5 = time_dist("planet1.rk5.dat")
time1, dist1 = time_dist("../planet1.dat")

# Compute the error on the planet distance
error_d0=np.max(np.abs((dist0-dist0_rk5)/dist0_rk5))
error_d1=np.max(np.abs((dist1-dist1_rk5)/dist1_rk5))


if plot:
    import matplotlib.pyplot as plt
    plt.close("all")
    fig, ax = plt.subplots(figsize=(6,6))
    ax.plot(time0_rk5, dist0_rk5, ls="-", c="k", label="dist 0 rk5")
    ax.plot(time0, dist0, ls="--", c="r", label="dist 0 new")
    ax.plot(time1_rk5, dist1_rk5, ls="-", c="b", label="dist 1 rk5")
    ax.plot(time1, dist1, ls="--", c="g", label="dist 1 new")
    ax.set(xlabel="time")
    fig.tight_layout()
    ax.legend(frameon=False)
    fig1, ax1 = plt.subplots(figsize=(6,6))
    ax1.plot(time0_active, active0, ls="-", c="k", label="isActive 0")
    ax1.plot(time1_active, active1, ls="--", c="r", label="isActive 1")
    ax1.axvline(x=time0_active[i0_active], c="k", ls=":")
    ax1.axvline(x=time1_active[i1_active], c="r", ls=":")
    ax1.set(xlabel="time", ylabel="planet activation")
    fig1.tight_layout()
    ax1.legend(frameon=False)
    plt.show()

print("Error_dist0=%e"%error_d0)
print("Error_dist1=%e"%error_d1)
print("active0[ip0]=%e"%active0[i0_active])
print("active1[ip1]=%e, active1[ip1-1]=%e"%(active1[i1_active],active1[i1_active-1]))
if error_d0<5.0e-10 and error_d1<5.0e-10 and active0[i0_active]==1.0 and active1[i1_active]==1.0 and active1[i1_active-1]==0.0:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
