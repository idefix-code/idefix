######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/HD/sod/pydefix_example.py
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
import pydefix as pdfx


def output(data):

    plt.figure()
    plt.plot(data.x[pdfx.IDIR], data.Vc[pdfx.VX1, 0, 0, :], label="VX1")
    plt.plot(data.x[pdfx.IDIR], data.Vc[pdfx.RHO, 0, 0, :], label="RHO")
    plt.legend()
    plt.show()
    # data.Vc[0,0,0,10] = 2.0


def initflow(data):

    # Initialize the flow
    data.Vc[pdfx.RHO, 0, 0, :] = 1.0
    data.Vc[pdfx.PRS, 0, 0, :] = 2.0
    data.Vc[pdfx.VX1, 0, 0, :] = np.sin(2.0 * np.pi * data.x[0][:])
