#!/usr/bin/env python
# coding: utf-8

import sys
import os
sys.path.append(os.getenv("IDEFIX_DIR"))
from pytools.vtk_io import readVTK

import numpy as np
import matplotlib.pyplot as plt

on = 3
ds = readVTK(f"../data.{on:04d}.vtk", geometry="polar")

Rmed = ds.x
phimed = (ds.y/(np.pi))%2*np.pi-np.pi
dR = np.ediff1d(ds.x)[0]

h0 = 0.05
flaring = 0.0

# Planet in a fixed circular orbit
Rp = 1.0
hp = h0*pow(Rp,flaring)

# See Bae+Zhu (2018, I)
# Dominant azimuthal wavenumber when
# perturbation driven by point mass perturber
mdom = int((1/2)*pow(hp,-1))

# Lindblad resonance radii
Rmp = pow(1+1/mdom,2/3)*Rp
Rmm = pow(1-1/mdom,2/3)*Rp

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def omega(
    R,
    *,
    keplerian=True,
    info=False,
):
    if info:
        print(f"keplerian = {keplerian}")
    if keplerian:
        return(pow(R,-3/2))

def soundSpeed(
    R,
    *,
    h0=h0,
    flaring=flaring,
    isotherm=True,
    info=False,
):
    if info:
        print(f"isotherm = {isotherm}")
        print(f"h0 = {h0:.2f}")
        print(f"flaring = {flaring:.2f}")
    if isotherm:
        return(h0*pow(R,flaring-1/2))

def _integral(
    RR,
    *,
    resonance=-1,
    m=mdom,
    Rp=Rp,
    h0=h0,
    flaring=flaring,
    info=False,
):
    if info:
        print(f"planet in Rp = {Rp:.2f}")
        print(f"dominant phi-wavenumber m = {m}")
    Rm = pow(1+np.sign(resonance)*1/m,2/3)*Rp
    data =     (omega(Rmed, info=info)/soundSpeed(Rmed, h0=h0, flaring=flaring, info=info))*pow(abs(pow(1-pow(Rmed/Rp,3/2),2)-1/m/m),1/2)
    k = find_nearest(Rmed, Rm)
    kk = find_nearest(Rmed, RR)
    if np.sign(resonance) > 0:
        integ = (
            np.nansum(
                (data * dR)[k : kk + 1],
                dtype="float64",
            )
        )
    else:
        integ = -(
            np.nansum(
                (data * dR)[kk : k + 1],
                dtype="float64",
            )
        )
    return integ

# Phase equation for m,n
#(phi-wavenumber,order of the spiral)
def phaseEquation(
    RR,
    *,
    resonance=-1,
    m=mdom,
    n=0,
    Rp=Rp,
    info=False,
):
    phieq = -np.sign(RR-Rp)*np.pi/4/m + 2*np.pi*n/m -_integral(
                RR,
                resonance=resonance,
                m=mdom,
                Rp=Rp,
                h0=h0,
                flaring=flaring,
                info=info,
            )
    return phieq

plot = False

RwkzMin = find_nearest(Rmed, 0.7)
RwkzMax = find_nearest(Rmed, 1.35)

spiralR = []
spiralP = []
spiralP_theo = []
for ir in range(RwkzMin, RwkzMax+1):
    if (Rmed[ir] < Rmm):
        rho = ds.data["RHO"][ir,:,0]
        spiralR.append(Rmed[ir])
        spiralP.append(phimed[find_nearest(rho,rho.max())])
        spiralP_theo.append(phaseEquation(Rmed[ir],resonance=-1))
    if (Rmed[ir] > Rmp):
        rho = ds.data["RHO"][ir,:,0]
        spiralR.append(Rmed[ir])
        spiralP.append(phimed[find_nearest(rho,rho.max())])
        spiralP_theo.append(phaseEquation(Rmed[ir],resonance=+1))

spiralR_theo = spiralR
spiralP = np.array(spiralP)
spiralP_theo = (np.array(spiralP_theo)/(np.pi))%2*np.pi-np.pi

if plot:
    fig, ax = plt.subplots(figsize=(5,4))
    ax.axvline(x=1.0, ls="--", c="k")
    ax.scatter(spiralR_theo,spiralP_theo, c="r",marker="+", label="theoretical")
    ax.scatter(spiralR,spiralP,c="k",marker="x", label=r"max $\rho$ simulation")
    fig.tight_layout()

    ax.legend(frameon=False)

    fig2, ax2 = plt.subplots()
    ax2.scatter(spiralR, 100*(spiralP-spiralP_theo)/spiralP_theo)
    fig2.tight_layout()

    plt.show()

mean = np.mean(100*(spiralP-spiralP_theo)/spiralP_theo)
std = np.std(100*(spiralP-spiralP_theo)/spiralP_theo)
error_mean = mean
error_dispersion = np.max([abs(mean+std),abs(mean-std)])

print("Mean Error=%.2f"%error_mean)
print("max(mean +/- std)=%.2f"%error_dispersion)
if error_mean<1.0 and error_dispersion<5.0:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
