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
planet_ref = datafile("planet0.ref.dat")
planet = datafile("../planet0.dat")

migration_ref = np.sqrt(planet_ref[1]**2+planet_ref[2]**2+planet_ref[3]**2)
migration = np.sqrt(planet[1]**2+planet[2]**2+planet[3]**2)

# Compute the error on the planet distance
error_d=np.max(np.abs((migration-migration_ref)/migration_ref))

torque_ref = datafile("tqwk0.ref.dat")
torque = datafile("../tqwk0.dat")

tq_tot_ref = torque_ref[1]+torque_ref[2]
tq_tot = torque[1]+torque[2]

# Compute the error on the planet torque
error_t=np.max(np.abs((tq_tot-tq_tot_ref)/tq_tot_ref))

V = readVTK("../data.0010.vtk", geometry="polar")
U = readVTK("data.0010.ref.vtk", geometry="polar")

# Compute the error on RHO
error_rho = np.mean(np.abs((V.data['RHO']-U.data['RHO'])/U.data['RHO']),axis=(0,1))

if plot:
    import matplotlib.pyplot as plt
    plt.close("all")
    fig, ax = plt.subplots(figsize=(6,6))
    ax.plot(planet_ref[-1], migration_ref, ls="-", c="k", label="dist ref")
    ax.plot(planet[-1], migration, ls="--", c="r", label="dist new")
    ax.set(xlabel="time")
    fig.tight_layout()
    ax.legend(frameon=False)
    fig2, ax2 = plt.subplots(figsize=(6,6))
    ax2.plot(torque_ref[-1], tq_tot_ref, ls="-", c="k", label="torque ref")
    ax2.plot(torque[-1], tq_tot, ls="--", c="r", label="torque new")
    ax2.set(xlabel="time")
    fig2.tight_layout()
    ax2.legend(frameon=False)
    fig3, ax3 = plt.subplots(ncols=2, figsize=(12,6))
    RR, PP = np.meshgrid(U.x, U.y%(2*np.pi)-np.pi, indexing="ij")
    imref = ax3[0].pcolormesh(RR, PP, U.data["VX1"][:,:,0], cmap="RdBu_r")
    ax3[0].set_title("ref")
    cbarref = fig3.colorbar(imref, ax=ax3[0])
    cbarref.set_label("VX1")
    RR2, PP2 = np.meshgrid(V.x, V.y%(2*np.pi)-np.pi, indexing="ij")
    imnew = ax3[1].pcolormesh(RR2, PP2, V.data["VX1"][:,:,0], cmap="RdBu_r")
    ax3[1].set_title("new")
    cbarnew = fig3.colorbar(imnew, ax=ax3[1])
    cbarnew.set_label("VX1")
    plt.show()

diff_torque = np.max(np.abs(tq_tot-tq_tot_ref))
print("diff_torque=%e"%diff_torque)
print("Error_dist=%e"%error_d)
# print("Error_torque=%e"%error_t)
print("Error_rho=%e"%error_rho)
# if error_d<5.0e-2 and error_t<5.0e-2 and error_rho<5.0e-2:
if error_d<5.0e-2 and diff_torque<1.0e-15 and error_rho<5.0e-2:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
