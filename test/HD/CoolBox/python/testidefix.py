#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Validate the tabulated optically-thin cooling module.

The CoolBox setup evolves a uniform box that cools from Tini down to the
temperature floor. The analysis file stores, for every snapshot:
    col0 : time
    col1 : mass-averaged temperature [K]
    col2 : instantaneous cooling rate  -de/dt / nH^2  [erg cm^3 s^-1]

We compare the measured cooling rate (col2) as a function of temperature
(col1) against the tabulated cooling curve used as input (cooltable.dat,
row 0 = temperature, row 1 = cooling rate), over the temperature range
[1.05e4, 1.95e6] K.

@author: dutta-alankar
"""

import sys
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-noplot",
                    default=False,
                    help="disable plotting",
                    action="store_true")
args, unknown = parser.parse_known_args()

# temperature range over which we validate the cooling rate
Tlo = 1.05e4
Thi = 1.95e6
tolerance = 5.0e-2

# tabulated cooling curve used as input by Idefix
table = np.loadtxt('../cooltable.dat')
Ttab = table[0, :]
Ltab = table[1, :]

# measured evolution
data = np.loadtxt('../analysis-time-temperature.dat')
Tsim = data[:, 1]
Lsim = data[:, 2]

# keep points with a meaningful (positive) cooling rate within the range
mask = (Tsim >= Tlo) & (Tsim <= Thi) & (Lsim > 0.0)
Tsim = Tsim[mask]
Lsim = Lsim[mask]

# reference cooling rate from the table at the measured temperatures (log-log)
Lref = np.power(10.0, np.interp(np.log10(Tsim), np.log10(Ttab), np.log10(Ltab)))

error = np.sqrt(np.mean((Lsim / Lref - 1.0) ** 2))

if not args.noplot:
    import matplotlib.pyplot as plt
    plt.close('all')
    plt.figure()
    plt.loglog(Ttab, Ltab, label='cooltable.dat')
    plt.loglog(Tsim, Lsim, '+', label='CoolBox')
    plt.xlabel('Temperature [K]')
    plt.ylabel(r'$\Lambda(T)$ [erg cm$^3$ s$^{-1}$]')
    plt.legend()
    print("Error=%e" % error)
    plt.ioff()
    plt.show()

if error < tolerance:
    print("SUCCESS")
    sys.exit(0)
else:
    print("Failed with error=%e" % error)
    sys.exit(1)
