#!/usr/bin/env python3
######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/utils/dumpImage/testme.py
#
# Last modified : 07/2026
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG / CNRS - 2023 - 2024)
# - Clément Robert <cr52@protonmail.com> (2026)
# - and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################
# Other contributors :
# - Gaylor Wafflard <gaylor.wafflard@univ-grenoble-alpes.fr> (IPAG / CNRS - 2024)
######################################################################################

"""

@author: glesur
"""

import os
import sys

sys.path.append(os.getenv("IDEFIX_DIR"))
import pytools.idfx_test as tst

test = tst.idfxTest(__file__)

test.configure()
test.compile()
# this test succeeds if it runs successfully
test.run()

test.mpi = True
test.configure()
test.compile()
# this test succeeds if it runs successfully
test.run()
