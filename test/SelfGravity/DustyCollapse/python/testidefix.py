import numpy as np
import os
import sys

sys.path.append(os.getenv("IDEFIX_DIR"))
from pytools.vtk_io import readVTK

G=6.6743e-11
M=1e16                 #mass of spheric clump
r0=1e4                 #radius of spheric clump
rmax = 3e4             #right edge of simulation box
d=M/(4*np.pi*r0**3/3)  #density of spheric clump
phi0 = -22.247666667   #gravitationnal potential at rmax

def potential(r):
    retv = np.empty_like(r)
    interior = (r<r0)
    retv[interior] = -4*np.pi*G*d*(3*r0**2-r[interior]**2)/6
    retv[~interior] = -4*np.pi*G*d*r0**3/(3*r[~interior])
    return retv


def test_potential():

    V = readVTK("../data.0001.vtk")

    r = V.r
    phi = np.squeeze(V.data["phiP"])
    pot = potential(r)

    np.testing.assert_allclose(pot, phi+phi0, rtol=2e-3)


if __name__ == "__main__":

    try:
        test_potential()
    except Exception:
        print("FAILURE!")
        sys.exit(1)
    else:
        print("SUCCESS!")
        sys.exit(0)
