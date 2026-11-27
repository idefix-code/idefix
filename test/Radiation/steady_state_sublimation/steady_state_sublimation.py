######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/Radiation/steady_state_sublimation/steady_state_sublimation.py
#
# Last modified : 09/2026
#
# Copyright(C) by :
# - Nicolas Scepi <nicolas.scepi@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
# - Sebastien Valat <sebastien.valat@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

import os

import numpy as np
from dump_io import readDump
from pydefix import *


def solve_hydrostatic(dr, dtheta, r, theta, vphi, T, Mdot):

    C_kb = 1.380649e-16
    c_mu = 1.6605390666e-24
    C_G = 6.674299999999999e-8
    Msol = 1.988409870698051e33
    GM = C_G * Msol
    mu = 2.35

    mid = int(theta.shape[1] / 2)
    nx1 = int(theta.shape[0])
    nx2 = int(theta.shape[1])
    nx3 = int(theta.shape[2])

    rho = np.zeros((nx1, nx2, nx3))
    P = np.zeros((nx1, nx2, nx3))

    R = r * np.sin(theta)

    Omega = vphi / R

    alphaMRI = 1.0e-1
    alphaDZ = 1.0e-3
    T_MRI = 900.0
    alpha_width = 25.0

    dens_floor = 1.0e-15 * np.ones(nx1)

    alpha = (alphaMRI - alphaDZ) * 0.5 * (
        1.0 - np.tanh((T_MRI - T) / alpha_width)
    ) + alphaDZ
    cs = np.sqrt(C_kb * T / (mu * c_mu))
    h = cs / Omega
    nu = alpha * cs**2 / Omega
    Sigma = Mdot / (3.0 * np.pi * nu[:, mid, :])
    rho[:, mid, :] = Sigma / (np.sqrt(2.0 * np.pi) * h[:, mid, :])

    P[:, mid, :] = rho[:, mid, :] * C_kb * T[:, mid, :] / (mu * c_mu)
    vphi[:, mid, 0] = np.sqrt(
        r[:, mid, 0]
        * np.diff(P[:, mid, 0], axis=0, append=2.0 * P[-1, mid, 0] - P[-2, mid, 0])
        / dr[:, mid]
        / rho[:, mid, 0]
        + GM / r[:, mid, 0]
    )

    for j in range(mid + 1, nx2):
        P[:, j, 0] = (
            P[:, j - 1, 0]
            + dtheta[:, j - 1]
            / np.tan(theta[:, j - 1, 0])
            * rho[:, j - 1, 0]
            * vphi[:, j - 1, 0] ** 2
        )
        rho[:, j, 0] = np.maximum(
            (P[:, j, 0] * mu * c_mu / (C_kb * T[:, j, 0])), dens_floor
        )
        P[:, j, 0] = rho[:, j, 0] * C_kb * T[:, j, 0] / (c_mu * mu)
        vphi[:, j, 0] = np.sqrt(
            r[:, j, 0]
            * np.diff(P[:, j, 0], append=2.0 * P[-1, j, 0] - P[-2, j, 0])
            / dr[:, j]
            / rho[:, j, 0]
            + GM / r[:, j, 0]
        )

    rho[:, :mid, 0] = rho[:, mid:, 0][:, ::-1]
    vphi[:, :mid, 0] = vphi[:, mid:, 0][:, ::-1]
    P[:, :mid, 0] = P[:, mid:, 0][:, ::-1]

    return rho, P, vphi


def init(data):
    dir_path = os.path.dirname(os.path.realpath(__file__))

    dump = readDump(dir_path + "/dump.ini.dmp")

    [r_grid, theta_grid, phi_grid] = np.meshgrid(
        dump.x1, dump.x2, dump.x3, indexing="ij"
    )
    dx1 = np.diff(dump.x1l, append=dump.x1r[-1])
    dx2 = np.diff(dump.x2l, append=dump.x2r[-1])
    dr_grid, dtheta_grid = np.meshgrid(dx1, dx2, indexing="ij")

    Msol = 1.988409870698051e33

    unit_velocity = 1.496e10
    unit_length = 1.496e13
    unit_density = 4.0e-28
    unit_energy = unit_density * unit_velocity**2
    KELVIN = 2.69164e12
    Msolyr_to_cgs = Msol / (3600.0 * 24.0 * 365.0)
    mu = 2.35

    Mdot = 1.0e-8 * Msolyr_to_cgs  ### in g

    T = (dump.data["Vc-PRS"] / dump.data["Vc-RHO"]) * KELVIN * mu

    vphi = dump.data["Vc-VX3"] * unit_velocity

    rho, P, vphi = solve_hydrostatic(
        dr_grid * unit_length,
        dtheta_grid,
        r_grid * unit_length,
        theta_grid,
        vphi,
        T,
        Mdot,
    )

    gbeg1 = data.gbeg[1] - data.nghost[1]
    gend1 = data.gend[1] - data.nghost[1]
    gbeg0 = data.gbeg[0] - data.nghost[0]
    gend0 = data.gend[0] - data.nghost[0]

    # Initialize the flow
    data.Vc[
        RHO, :, data.nghost[1] : -data.nghost[1], data.nghost[0] : -data.nghost[0]
    ] = rho.T[:, gbeg1:gend1, gbeg0:gend0] / unit_density
    data.Vc[
        PRS, :, data.nghost[1] : -data.nghost[1], data.nghost[0] : -data.nghost[0]
    ] = P.T[:, gbeg1:gend1, gbeg0:gend0] / unit_energy
    data.Vc[
        VX1, :, data.nghost[1] : -data.nghost[1], data.nghost[0] : -data.nghost[0]
    ] = 0.0
    data.Vc[
        VX2, :, data.nghost[1] : -data.nghost[1], data.nghost[0] : -data.nghost[0]
    ] = 0.0
    data.Vc[
        VX3, :, data.nghost[1] : -data.nghost[1], data.nghost[0] : -data.nghost[0]
    ] = vphi.T[:, gbeg1:gend1, gbeg0:gend0] / unit_velocity
