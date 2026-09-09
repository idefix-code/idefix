######################################################################################
# Idefix MHD astrophysical code
#
# Source file pytools/tests/test_idfx_test_run.py
#
# Last modified : 07/2026
#
# Copyright(C) by :
# - Clément Robert <cr52@protonmail.com> (2026)
# - Sébastien Valat <sebastien.valat@univ-grenoble-alpes.fr> (IPAG / CNRS - 2026)
# - and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

import os

import pytest

from ..idfx_test_run import IdexPytestRunner


def test_genTests():
    # build runner
    runner = IdexPytestRunner(__file__)

    # dir
    dir = os.path.dirname(__file__)

    # generate
    result = runner.genTests()
    assert result == [
        pytest.param(
            {
                "dumpname": "dump.0001.dmp",
                "noplot": True,
                "reconstruction": 2,
                "tolerance": 1e-14,
                "ini": "idefix.ini",
                "testfile": dir + "/test/pb1/testme.json",
                "testname": "pb1",
            },
            marks=(),
            id="pb1-idefix.ini",
        ),
        pytest.param(
            {
                "dumpname": "dump.0001.dmp",
                "noplot": True,
                "reconstruction": 2,
                "tolerance": 1e-14,
                "ini": "idefix-implicit.ini",
                "testfile": dir + "/test/pb1/testme.json",
                "testname": "pb1",
            },
            marks=(),
            id="pb1-idefix-implicit.ini",
        ),
        pytest.param(
            {
                "dumpname": "dump.0001.dmp",
                "noplot": True,
                "reconstruction": 2,
                "tolerance": 1e-14,
                "ini": "idefix.ini",
                "testfile": dir + "/test/pb2/testme.json",
                "testname": "pb2",
            },
            marks=(),
            id="pb2-idefix.ini-noplot",
        ),
        pytest.param(
            {
                "dumpname": "dump.0001.dmp",
                "noplot": True,
                "reconstruction": 2,
                "tolerance": 1e-14,
                "ini": "idefix-implicit.ini",
                "testfile": dir + "/test/pb2/testme.json",
                "testname": "pb2",
            },
            marks=(),
            id="pb2-idefix-implicit.ini-noplot",
        ),
        pytest.param(
            {
                "dumpname": "dump.0001.dmp",
                "noplot": False,
                "reconstruction": 2,
                "tolerance": 1e-14,
                "ini": "idefix.ini",
                "testfile": dir + "/test/pb2/testme.json",
                "testname": "pb2",
            },
            marks=(),
            id="pb2-idefix.ini",
        ),
        pytest.param(
            {
                "dumpname": "dump.0001.dmp",
                "noplot": False,
                "reconstruction": 2,
                "tolerance": 1e-14,
                "ini": "idefix-implicit.ini",
                "testfile": dir + "/test/pb2/testme.json",
                "testname": "pb2",
            },
            marks=(),
            id="pb2-idefix-implicit.ini",
        ),
    ]
