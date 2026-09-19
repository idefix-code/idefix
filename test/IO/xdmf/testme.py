#!/usr/bin/env python3
######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/IO/xdmf/testme.py
#
# Last modified : 07/2026
#
# Copyright(C) by :
# - Clément Robert <clement.robert@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

"""

@author: glesur
"""

import os
import sys

sys.path.append(os.getenv("IDEFIX_DIR"))

import pytools.idfx_test as tst

name = "dump.0001.dmp"

test = tst.idfxTest(__file__)


def check_xdmf_exists():
    # verify that the expected XDMF sidecar and data file exist.
    xdmf_file = "data.0001.flt.xmf"
    h5_file = "data.0001.flt.h5"

    if not os.path.exists(xdmf_file):
        print("Missing expected XDMF output file: {}".format(xdmf_file))
        sys.exit(1)
    if not os.path.exists(h5_file):
        print("Missing expected XDMF data file: {}".format(h5_file))
        sys.exit(1)


if not test.dec:
    test.dec = ["2", "2", "2"]

if not test.check:
    test.vectPot = False
    test.single = False
    test.reconstruction = 2
    test.mpi = False
    # Only check that the test runs.
    test.configure(definitionFile="definitions.hpp")
    test.compile()
    test.run(inputFile="idefix.ini")

check_xdmf_exists()
