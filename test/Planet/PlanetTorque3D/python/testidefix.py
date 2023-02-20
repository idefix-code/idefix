#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  21 12:04:41 2021

@author: gwafflard
"""

import os
import sys
import re
sys.path.append(os.getenv("IDEFIX_DIR"))
from pytools.vtk_io import readVTK
import numpy as np

def datafile(filename, *, directory=""):
    fullpath = os.path.join(directory, filename)
    with open(fullpath) as f1:
        data = f1.readlines()
    y = [[v for v in re.split(r"[\t ]+", r)] for r in data]
    columns = np.array(y, dtype="float64").T
    return(columns)

plot = False

torque_ref = datafile("tqwk0.ref.dat")
torque = datafile("../tqwk0.dat")

tq_tot_ref = torque_ref[1]+torque_ref[2]
tq_tot = torque[1]+torque[2]
tqh_tot_ref = torque_ref[3]+torque_ref[4]
tqh_tot = torque[3]+torque[4]

# Compute the error on the planet torque
error_t=np.max(np.abs((tq_tot-tq_tot_ref)/tq_tot_ref))
error_th=np.max(np.abs((tqh_tot-tqh_tot_ref)/tqh_tot_ref))

V = readVTK("../data.0009.vtk", geometry="spherical")
U = readVTK("data.0009.ref.vtk", geometry="spherical")

# Compute the error on RHO
error_rho = np.mean(np.abs((V.data['RHO']-U.data['RHO'])/U.data['RHO']),axis=(0,1,2))

if plot:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6,6))
    ax.plot(torque_ref[-1], tq_tot_ref, ls="-", c="k", label="torque ref")
    ax.plot(torque[-1], tq_tot, ls="--", c="r", label="torque new")
    ax.plot(torque_ref[-1], tqh_tot_ref, ls="-", c="b", label="torque nohill ref")
    ax.plot(torque[-1], tqh_tot, ls="--", c="g", label="torque nohill new")
    ax.set(xlabel="time")
    fig.tight_layout()
    ax.legend(frameon=False)
    fig3, ax3 = plt.subplots(ncols=2, figsize=(12,6))
    RR, PP = np.meshgrid(U.r, U.phi%(2*np.pi)-np.pi, indexing="ij")
    imref = ax3[0].pcolormesh(RR, PP, U.data["RHO"][:,len(U.theta)//2,:], cmap="viridis")
    ax3[0].set_title("ref")
    cbarref = fig3.colorbar(imref, ax=ax3[0])
    cbarref.set_label("RHO")
    RR2, PP2 = np.meshgrid(V.r, V.phi%(2*np.pi)-np.pi, indexing="ij")
    imnew = ax3[1].pcolormesh(RR2, PP2, V.data["RHO"][:,len(U.theta)//2,:], cmap="viridis")
    ax3[1].set_title("new")
    cbarnew = fig3.colorbar(imnew, ax=ax3[1])
    cbarnew.set_label("RHO")
    plt.show()

diff_torque = np.max(np.abs(tq_tot-tq_tot_ref))
diff_torqueh = np.max(np.abs(tqh_tot-tqh_tot_ref))
print("diff_torque=%e"%diff_torque)
print("diff_torque_nohill=%e"%diff_torqueh)
# print("Error_torque=%e"%error_t)
# print("Error_torque_nohill=%e"%error_th)
print("Error_rho=%e"%error_rho)
# if error_t<5.0e-2 and error_th<5.0e-2 and error_rho<5.0e-2:
if diff_torque<1.0e-15 and diff_torqueh<1.0e-15 and error_rho<5.0e-2:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
