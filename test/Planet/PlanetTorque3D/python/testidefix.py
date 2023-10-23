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

def datafile(filename):
    return np.loadtxt(filename, dtype="float64").T

plot = False

torque_ref = datafile("tqwk0.ref.dat")
torque = datafile("../tqwk0.dat")

# Factor 2 here because the reference has been computed for half the disk
# (hence the torque was half)
tq_tot_ref = 2*(torque_ref[1]+torque_ref[2])
tq_tot = torque[1]+torque[2]
tqh_tot_ref = 2*(torque_ref[3]+torque_ref[4])
tqh_tot = torque[3]+torque[4]

# Compute the error on the planet torque
error_t=np.max(np.abs((tq_tot-tq_tot_ref)/tq_tot_ref))
error_th=np.max(np.abs((tqh_tot-tqh_tot_ref)/tqh_tot_ref))

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

    plt.show()

diff_torque = np.max(np.abs(tq_tot-tq_tot_ref))
diff_torqueh = np.max(np.abs(tqh_tot-tqh_tot_ref))
print("diff_torque=%e"%diff_torque)
print("diff_torque_nohill=%e"%diff_torqueh)
# print("Error_torque=%e"%error_t)
# print("Error_torque_nohill=%e"%error_th)
# if error_t<5.0e-2 and error_th<5.0e-2 and error_rho<5.0e-2:
if diff_torque<1.0e-12 and diff_torqueh<1.0e-13:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
