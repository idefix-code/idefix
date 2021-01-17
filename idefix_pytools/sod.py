"""
Created on Thu Mar  5 11:27:16 2020

@author: glesur
Source from
https://github.com/ibackus/sod-shocktube/blob/master/sod.py
"""

import numpy as np
import scipy
import scipy.optimize


def sound_speed(gamma, pressure, density, dustFrac=0.):
    """
    Calculate sound speed, scaled by the dust fraction according to:

        .. math::
            \widetilde{c}_s = c_s \sqrt{1 - \epsilon}

    Where :math:`\epsilon` is the dustFrac
    """
    scale = np.sqrt(1 - dustFrac)
    return np.sqrt(gamma * pressure/ density) * scale

def shock_tube_function(p4, p1, p5, rho1, rho5, gamma, dustFrac=0.):
    """
    Shock tube equation
    """
    z = (p4 / p5 - 1.)
    c1 = sound_speed(gamma, p1, rho1, dustFrac)
    c5 = sound_speed(gamma, p5, rho5, dustFrac)

    gm1 = gamma - 1.
    gp1 = gamma + 1.
    g2 = 2. * gamma

    fact = gm1 / g2 * (c5 / c1) * z / np.sqrt(1. + gp1 / g2 * z)
    fact = (1. - fact) ** (g2 / gm1)

    return p1 * fact - p4


def calculate_regions(pl, ul, rhol, pr, ur, rhor, gamma=1.4, dustFrac=0.):
    """
    Compute regions
    :rtype : tuple
    :return: returns p, rho and u for regions 1,3,4,5 as well as the shock speed
    """
    # if pl > pr...
    rho1 = rhol
    p1 = pl
    u1 = ul
    rho5 = rhor
    p5 = pr
    u5 = ur

    # unless...
    if pl < pr:
        rho1 = rhor
        p1 = pr
        u1 = ur
        rho5 = rhol
        p5 = pl
        u5 = ul

    # solve for post-shock pressure
    p4 = scipy.optimize.fsolve(shock_tube_function, p1, (p1, p5, rho1, rho5, gamma))[0]

    # compute post-shock density and velocity
    z = (p4 / p5 - 1.)
    c5 = sound_speed(gamma, p5, rho5, dustFrac)

    gm1 = gamma - 1.
    gp1 = gamma + 1.
    gmfac1 = 0.5 * gm1 / gamma
    gmfac2 = 0.5 * gp1 / gamma

    fact = np.sqrt(1. + gmfac2 * z)

    u4 = c5 * z / (gamma * fact)
    rho4 = rho5 * (1. + gmfac2 * z) / (1. + gmfac1 * z)

    # shock speed
    w = c5 * fact

    # compute values at foot of rarefaction
    p3 = p4
    u3 = u4
    rho3 = rho1 * (p3 / p1)**(1. / gamma)
    return (p1, rho1, u1), (p3, rho3, u3), (p4, rho4, u4), (p5, rho5, u5), w


def calc_positions(pl, pr, region1, region3, w, xi, t, gamma, dustFrac=0.):
    """
    :return: tuple of positions in the following order ->
            Head of Rarefaction: xhd,  Foot of Rarefaction: xft,
            Contact Discontinuity: xcd, Shock: xsh
    """
    p1, rho1 = region1[:2]  # don't need velocity
    p3, rho3, u3 = region3
    c1 = sound_speed(gamma, p1, rho1, dustFrac)
    c3 = sound_speed(gamma, p3, rho3, dustFrac)

    if pl > pr:
        xsh = xi + w * t
        xcd = xi + u3 * t
        xft = xi + (u3 - c3) * t
        xhd = xi - c1 * t
    else:
        # pr > pl
        xsh = xi - w * t
        xcd = xi - u3 * t
        xft = xi - (u3 - c3) * t
        xhd = xi + c1 * t

    return xhd, xft, xcd, xsh


def region_states(pl, pr, region1, region3, region4, region5):
    """
    :return: dictionary (region no.: p, rho, u), except for rarefaction region
    where the value is a string, obviously
    """
    if pl > pr:
        return {'Region 1': region1,
                'Region 2': 'RAREFACTION',
                'Region 3': region3,
                'Region 4': region4,
                'Region 5': region5}
    else:
        return {'Region 1': region5,
                'Region 2': region4,
                'Region 3': region3,
                'Region 4': 'RAREFACTION',
                'Region 5': region1}


def create_arrays(pl, pr, xl, xr, positions, state1, state3, state4, state5,
                  npts, gamma, t, xi, dustFrac=0.):
    """
    :return: tuple of x, p, rho and u values across the domain of interest
    """
    xhd, xft, xcd, xsh = positions
    p1, rho1, u1 = state1
    p3, rho3, u3 = state3
    p4, rho4, u4 = state4
    p5, rho5, u5 = state5
    gm1 = gamma - 1.
    gp1 = gamma + 1.

    x_arr = np.linspace(xl, xr, npts)
    rho = np.zeros(npts, dtype=float)
    p = np.zeros(npts, dtype=float)
    u = np.zeros(npts, dtype=float)
    c1 = sound_speed(gamma, p1, rho1, dustFrac)
    if pl > pr:
        for i, x in enumerate(x_arr):
            if x < xhd:
                rho[i] = rho1
                p[i] = p1
                u[i] = u1
            elif x < xft:
                u[i] = 2. / gp1 * (c1 + (x - xi) / t)
                fact = 1. - 0.5 * gm1 * u[i] / c1
                rho[i] = rho1 * fact ** (2. / gm1)
                p[i] = p1 * fact ** (2. * gamma / gm1)
            elif x < xcd:
                rho[i] = rho3
                p[i] = p3
                u[i] = u3
            elif x < xsh:
                rho[i] = rho4
                p[i] = p4
                u[i] = u4
            else:
                rho[i] = rho5
                p[i] = p5
                u[i] = u5
    else:
        for i, x in enumerate(x_arr):
            if x < xsh:
                rho[i] = rho5
                p[i] = p5
                u[i] = -u1
            elif x < xcd:
                rho[i] = rho4
                p[i] = p4
                u[i] = -u4
            elif x < xft:
                rho[i] = rho3
                p[i] = p3
                u[i] = -u3
            elif x < xhd:
                u[i] = -2. / gp1 * (c1 + (xi - x) / t)
                fact = 1. + 0.5 * gm1 * u[i] / c1
                rho[i] = rho1 * fact ** (2. / gm1)
                p[i] = p1 * fact ** (2. * gamma / gm1)
            else:
                rho[i] = rho1
                p[i] = p1
                u[i] = -u1

    return x_arr, p, rho, u


def solve(left_state, right_state, geometry, t, gamma=1.4, npts=500,
          dustFrac=0.):
    """
    Solves the Sod shock tube problem (i.e. riemann problem) of discontinuity
    across an interface.

    Parameters
    ----------
    left_state, right_state: tuple
        A tuple of the state (pressure, density, velocity) on each side of the
        shocktube barrier for the ICs.  In the case of a dusty-gas, the density
        should be the gas density.
    geometry: tuple
        A tuple of positions for (left boundary, right boundary, barrier)
    t: float
        Time to calculate the solution at
    gamma: float
        Adiabatic index for the gas.
    npts: int
        number of points for array of pressure, density and velocity
    dustFrac: float
        Uniform fraction for the gas, between 0 and 1.

    Returns
    -------
    positions: dict
        Locations of the important places (rarefaction wave, shock, etc...)
    regions: dict
        constant pressure, density and velocity states in distinct regions
    values: dict
        Arrays of pressure, density, and velocity as a function of position.
        The density ('rho') is the gas density, which may differ from the
        total density in a dusty-gas.
        Also calculates the specific internal energy
    """

    pl, rhol, ul = left_state
    pr, rhor, ur = right_state
    xl, xr, xi = geometry

    # basic checking
    if xl >= xr:
        print('xl has to be less than xr!')
        exit()
    if xi >= xr or xi <= xl:
        print('xi has in between xl and xr!')
        exit()

    # calculate regions
    region1, region3, region4, region5, w = \
        calculate_regions(pl, ul, rhol, pr, ur, rhor, gamma, dustFrac)

    regions = region_states(pl, pr, region1, region3, region4, region5)

    # calculate positions
    x_positions = calc_positions(pl, pr, region1, region3, w, xi, t, gamma,
                                 dustFrac)

    pos_description = ('Head of Rarefaction', 'Foot of Rarefaction',
                       'Contact Discontinuity', 'Shock')
    positions = dict(zip(pos_description, x_positions))

    # create arrays
    x, p, rho, u = create_arrays(pl, pr, xl, xr, x_positions,
                                 region1, region3, region4, region5,
                                 npts, gamma, t, xi, dustFrac)

    energy = p/(rho * (gamma - 1.0))
    rho_total = rho/(1.0 - dustFrac)
    val_dict = {'x':x, 'p':p, 'rho':rho, 'u':u, 'energy':energy,
                'rho_total':rho_total}

    return positions, regions, val_dict
