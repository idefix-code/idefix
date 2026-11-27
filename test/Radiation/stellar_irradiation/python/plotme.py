#!/usr/bin/env python3
######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/Radiation/stellar_irradiation/python/plotme.py
#
# Last modified : 11/2026
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

import argparse
import os
import sys
import astropy.units as u
import astropy.constants as c

sys.path.append(os.getenv("IDEFIX_DIR"))

import matplotlib.pyplot as plt
import numpy as np
import inifix

from pytools import sod
from pytools.vtk_io import readVTK


def main():
    parser = argparse.ArgumentParser(description="Plot radiation shock data.")
    parser.add_argument("-i", "--ini", type=str,
                        help="Path to the .ini file opened by inifix.")
    args = parser.parse_args()

    ini = inifix.load(args.ini)

    V = readVTK("../data.0001.vtk")

    unit_density = ini["Units"]["density"] * u.g/u.cm**3
    unit_length = ini["Units"]["length"]*u.cm
    unit_velocity = ini["Units"]["velocity"]*u.cm/u.s
    unit_mu = ini["Hydro"]["mu"]
    unit_energy = unit_density * unit_velocity**2
    unit_kelvin = u.u * unit_velocity**2 / c.k_B
    unit_ar = 4 * c.sigma_sb / c.c

    Tgas = ((V.data["PRS"] * unit_energy) / (V.data["RHO"] * unit_density / unit_mu / u.u) / c.k_B).to(u.K)

    Rl = V.r[:,None]*np.sin(V.theta[None,:])
    Zl = V.r[:,None]*np.cos(V.theta[None,:])

    plt.figure(1)
    plt.pcolormesh(Rl, Zl/Rl, np.log10(Tgas[:,:,0].value), shading="auto",cmap="inferno")
    plt.xscale("log")
    plt.ylim([-1,1])
    plt.xlim([0.9,1000])

    plt.colorbar(label="log(Tgas) [K]")
    

    plt.figure(2)
    nmid=V.theta.size//2
    plt.loglog(V.r, Tgas[:,nmid,0].value, label="Tgas")
    plt.xlabel("r [AU]")
    plt.ylabel("Tgas [K]")

    plt.figure(3)
    nmid=np.argwhere(V.r>2)[0][0]
    plt.plot(V.theta-np.pi/2, Tgas[nmid,:,0].value, label="Tgas")
    plt.xlabel(r"$\theta-\pi/2$")
    plt.ylabel("Tgas [K]")


    plt.show()

if __name__ == "__main__":
    main()
