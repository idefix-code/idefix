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

def separation(resonance):
    a = int(resonance[-1])
    b = int(resonance[0])
    return pow(a/b,2/3)

def datafile(filename):
    return np.loadtxt(filename, dtype="float64").T

planet0 = datafile("../planet0.dat")
planet1 = datafile("../planet1.dat")

resonance = "3:2"
plot_orbit = False

theta1 = np.arctan2(planet1[2],planet1[1])
apl1 = 1.0
ecc1 = 0.0
rpl1 = (apl1*(1-ecc1**2)/(1-ecc1*np.cos(theta1)))
xpl1 = rpl1*np.cos(theta1)
ypl1 = rpl1*np.sin(theta1)

theta0 = np.arctan2(planet0[2],planet0[1])
apl0 = separation(resonance)
ecc0 = 0.2
rpl0 = (apl0*(1-ecc0**2)/(1+ecc0*np.cos(theta0)))
xpl0 = rpl0*np.cos(theta0)
ypl0 = rpl0*np.sin(theta0)

if plot_orbit:
    import matplotlib.pyplot as plt
    plt.close("all")
    fig, ax = plt.subplots(figsize=(6,6))
    ax.scatter(planet0[1], planet0[2], ls="-", c="k")
    ax.scatter(planet1[1], planet1[2], ls="--", c="r")
    ax.scatter(0,0,marker="+")
    ax.plot(xpl0, ypl0, c="b")
    ax.plot(xpl1, ypl1, c="g")
    ax.set(xlabel="x", ylabel="y", aspect="equal", xlim=(-1.2,1.2), ylim=(-1.2,1.2))
    # ax.set(xlabel="x", ylabel="y", aspect="equal", xlim=(-2.5,2.5), ylim=(-2.5,2.5))
    ax.set_title(r"m$_{\rm sat}=0$, e=%.1f, n/n$_s$=%s (a0=%.3f)"%(ecc0,resonance,separation(resonance=resonance)), pad=20)
    # plt.savefig(f"resonance{resonance}_e{ecc}.png", dpi=200)
    fig.tight_layout()

    fig2, ax2 = plt.subplots(figsize=(6,6))
    ax2.scatter(theta0, (rpl0-np.sqrt(planet0[1]**2 + planet0[2]**2))/(rpl0), c="b", s=200)
    ax2.scatter(theta1, (rpl1-np.sqrt(planet1[1]**2 + planet1[2]**2))/(rpl1), c="g", s=200)

    plt.show()

rpl0_sim = np.sqrt(planet0[1]**2 + planet0[2]**2)
rpl1_sim = np.sqrt(planet1[1]**2 + planet1[2]**2)

mean0 = np.mean(100*(rpl0-rpl0_sim)/rpl0)
mean1 = np.mean(100*(rpl1-rpl1_sim)/rpl1)
std0 = np.std(100*(rpl0-rpl0_sim)/rpl0)
std1 = np.std(100*(rpl1-rpl1_sim)/rpl1)
error_mean0 = mean0
error_mean1 = mean1
error_dispersion0 = np.max([abs(mean0+std0),abs(mean0-std0)])
error_dispersion1 = np.max([abs(mean1+std1),abs(mean1-std1)])

print("Mean Error (planet0)=%s"%error_mean0)
print("Mean Error (planet1)=%s"%error_mean1)
print("max(mean +/- std) (planet0)=%s"%error_dispersion0)
print("max(mean +/- std) (planet1)=%s"%error_dispersion1)
if error_mean0<1.0e-10 and error_dispersion0<1.0e-9 and error_mean1<1.0e-10 and error_dispersion1<1.0e-9:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
