#!/usr/bin/env python3

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

    V = readVTK("../data.0001.vtk", geometry="cartesian")

    unit_density = ini["Units"]["density"] * u.g/u.cm**3
    unit_length = ini["Units"]["length"]*u.cm
    unit_velocity = ini["Units"]["velocity"]*u.cm/u.s
    unit_mu = ini["Hydro"]["mu"]
    unit_energy = unit_density * unit_velocity**2
    unit_kelvin = u.u * unit_velocity**2 / c.k_B
    unit_ar = 4 * c.sigma_sb / c.c

    Tgas = ((V.data["PRS"] * unit_energy) / (V.data["RHO"] * unit_density / unit_mu / u.u) / c.k_B).to(u.K)
    Trad = (((V.data["Rad0_RHO"] * unit_energy) / (unit_ar))**(1/4)).to(u.K)

    Tplus = np.amax(Tgas[:, 0, 0].value)
    imax=np.argmax(Tgas[:, 0, 0].value)
    T2 = Tgas[0, 0, 0].value
    Tminus = Tgas[imax+1, 0, 0].value

    print(f"T+={Tplus:.2e} T-={Tminus:.2e} T_2={T2:.2e}")


    plt.figure(1)
    plt.plot(V.x, Tgas[:, 0, 0].value, label="Tgas",marker="o")
    plt.plot(V.x, Trad[:, 0, 0].value, label="Trad")
    plt.xlim([0,5])
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
