######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/IO/pydefix/pydefix_example.py
#
# Last modified : 07/2026
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2024)
# - Clément Robert <clement.robert@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

import matplotlib.pyplot as plt
import numpy as np
from pydefix import *


def test_input_parameters_access():
    """Test that input_parameters are accessible from Python."""
    # Access the global input_parameters dict
    params = input_parameters

    # Verify it's a dict
    assert isinstance(params, dict), "input_parameters should be a dict"

    # Verify we can access Grid parameters
    assert "Grid" in params, "Grid block should be in input_parameters"
    grid_params = params["Grid"]
    assert isinstance(grid_params, dict), "Grid should be a dict"

    # Verify we can access specific parameters
    assert "X1-grid" in grid_params, "X1-grid should be in Grid params"
    x1_grid = grid_params["X1-grid"]
    assert isinstance(x1_grid, list), "X1-grid should be a list"
    assert len(x1_grid) > 0, "X1-grid should have values"

    # Verify we can read actual values from the ini file
    # From idefix.ini: X1-grid is "1  0.0  256  u  1.0"
    assert x1_grid[0] == "1", f"Expected '1', got {x1_grid[0]}"
    assert x1_grid[1] == "0.0", f"Expected '0.0', got {x1_grid[1]}"
    assert x1_grid[2] == "256", f"Expected '256', got {x1_grid[2]}"

    print("✓ input_parameters access test passed!")


# The output function
# the only argument is dataBlockHost python object, wrapping a dataBlockHost Idefix object
def output(data, grid, n):
    # Note: if MPI is not enabled, GatherIdefixArray still works (but merely does a local copy)
    pressure = GatherIdefixArray(
        data.Vc[PRS, :, :, :], data, broadcast=False, keepBoundaries=True
    )
    # only process #0 performs the output
    if prank == 0:
        plt.close()
        plt.figure()
        plt.pcolormesh(
            grid.x[IDIR],
            grid.x[JDIR],
            pressure[0, :, :],
            label="PRS",
            vmin=0.02,
            vmax=0.5,
            cmap="plasma",
        )
        plt.title("t=%.2f" % data.t)
        plt.colorbar()
        plt.savefig("PRS.%.4d.png" % n)


def initflow(data):
    # Test input_parameters access on first call
    test_input_parameters_access()

    # Field amplitude
    B0 = 1 / np.sqrt(4 * np.pi)

    [z, y, x] = np.meshgrid(data.x[KDIR], data.x[JDIR], data.x[IDIR], indexing="ij")

    # Initialize the flow
    data.Vc[RHO, :, :, :] = 25 / (36 * np.pi)
    data.Vc[PRS, :, :, :] = 5 / (12 * np.pi)
    data.Vc[VX1, :, :, :] = -np.sin(2 * np.pi * y)
    data.Vc[VX2, :, :, :] = np.sin(2 * np.pi * x)

    data.Vs[BX1s, :, :-1, :-1] = -B0 * np.sin(2 * np.pi * y)
    data.Vs[BX2s, :, :-1, :-1] = B0 * np.sin(4 * np.pi * x)
