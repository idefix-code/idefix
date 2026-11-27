#!/usr/bin/env python3
######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/Radiation/pulse/testme.py
#
# Last modified : 11/2026
#
# Copyright(C) by :
# - Nicolas Scepi <nicolas.scepi@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
# - Sebastien Valat <sebastien.valat@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

"""

@author: glesur
"""

import os
import shutil
import sys

sys.path.append(os.getenv("IDEFIX_DIR"))

import pytools.idfx_test as tst

name = "dump.0001.dmp"

tolerance = 1e-10


def testMe(test):

    deffiles = [
        "definitions-sph1D.hpp",
        "definitions-sph2D.hpp",
        "definitions-sph3D.hpp",
        "definitions-cart3D.hpp",
    ]

    # loop on all the ini files for this test
    for definition in deffiles:
        inifiles = []
        for solver in ["hll", "hllc"]:
            for opt in ["thin", "thick"]:
                if "cart" in definition:
                    if "1D" in definition:
                        ini_name = (
                            "idefix-" + str(solver) + "-cart1D-" + str(opt) + ".ini"
                        )
                        shutil.copy(
                            "idefix-" + str(solver) + "-cart-" + str(opt) + ".ini",
                            str(ini_name),
                        )
                        inifiles += [str(ini_name)]
                        test.dec = ["2"]
                    elif "2D" in definition:
                        ini_name = (
                            "idefix-" + str(solver) + "-cart2D-" + str(opt) + ".ini"
                        )
                        shutil.copy(
                            "idefix-" + str(solver) + "-cart-" + str(opt) + ".ini",
                            str(ini_name),
                        )
                        inifiles += [str(ini_name)]
                        test.dec = ["2", "2"]
                    elif "3D" in definition:
                        ini_name = (
                            "idefix-" + str(solver) + "-cart3D-" + str(opt) + ".ini"
                        )
                        shutil.copy(
                            "idefix-" + str(solver) + "-cart-" + str(opt) + ".ini",
                            str(ini_name),
                        )
                        inifiles += [str(ini_name)]
                        test.dec = ["2", "2", "1"]
                elif "sph" in definition:
                    if "1D" in definition:
                        ini_name = (
                            "idefix-" + str(solver) + "-sph1D-" + str(opt) + ".ini"
                        )
                        shutil.copy(
                            "idefix-" + str(solver) + "-sph-" + str(opt) + ".ini",
                            str(ini_name),
                        )
                        inifiles += [str(ini_name)]
                        test.dec = ["2"]
                    elif "2D" in definition:
                        ini_name = (
                            "idefix-" + str(solver) + "-sph2D-" + str(opt) + ".ini"
                        )
                        shutil.copy(
                            "idefix-" + str(solver) + "-sph-" + str(opt) + ".ini",
                            str(ini_name),
                        )
                        inifiles += [str(ini_name)]
                        test.dec = ["2", "2"]
                    elif "3D" in definition:
                        if opt == "thick":
                            continue
                        ini_name = (
                            "idefix-" + str(solver) + "-sph3D-" + str(opt) + ".ini"
                        )
                        shutil.copy(
                            "idefix-" + str(solver) + "-sph-" + str(opt) + ".ini",
                            str(ini_name),
                        )
                        inifiles += [str(ini_name)]
                        test.dec = ["2", "2", "1"]

        test.configure(definitionFile=definition)
        test.compile()
        for ini in inifiles:
            test.run(inputFile=ini)
            if test.init and not test.mpi:
                test.makeReference(filename=name)
            test.standardTest()
            test.nonRegressionTest(filename=name, tolerance=tolerance)
            os.remove(ini)


test = tst.idfxTest(__file__)

if not test.all:
    if test.check:
        test.checkOnly(filename=name, tolerance=tolerance)
    else:
        testMe(test)
else:
    for rec in [2]:
        test.vectPot = False
        test.single = False
        test.reconstruction = rec
        test.mpi = False
        testMe(test)

    if "1D" not in test.definitions:
        test.mpi = True
        testMe(test)
