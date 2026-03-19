#####################################################################################
# Idefix MHD astrophysical code
# Copyright(C) Sébastien Valat <sebastien.valat@univ-grenoble-alpes.fr>
# and other code contributors
# Licensed under CeCILL 2.1 License, see COPYING for more information
#####################################################################################

from ..idfx_test_run import IdexPytestRunner
import pytest
import os

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
        'dumpname': 'dump.0001.dmp',
        'noplot': True,
        'reconstruction': 2,
        'tolerance': 1e-14,
        'ini': 'idefix.ini',
        'testfile': dir + '/test/pb1/testme.json',
        'testname': 'pb1'
      }, marks=(), id='pb1-idefix.ini'
    ),
    pytest.param(
      {
        'dumpname': 'dump.0001.dmp',
        'noplot': True,
        'reconstruction': 2,
        'tolerance': 1e-14,
        'ini': 'idefix-implicit.ini',
        'testfile': dir + '/test/pb1/testme.json',
        'testname': 'pb1'
      }, marks=(), id='pb1-idefix-implicit.ini'
    ),
    pytest.param(
      {
        'dumpname': 'dump.0001.dmp',
        'noplot': True,
        'reconstruction': 2,
        'tolerance': 1e-14,
        'ini': 'idefix.ini',
        'testfile': dir + '/test/pb2/testme.json',
        'testname': 'pb2'
      }, marks=(), id='pb2-idefix.ini-noplot'
    ),
    pytest.param(
      {
        'dumpname': 'dump.0001.dmp',
        'noplot': True,
        'reconstruction': 2,
        'tolerance': 1e-14,
        'ini': 'idefix-implicit.ini',
        'testfile': dir + '/test/pb2/testme.json',
        'testname': 'pb2'
      }, marks=(), id='pb2-idefix-implicit.ini-noplot'
    ),
     pytest.param(
      {
        'dumpname': 'dump.0001.dmp',
        'noplot': False,
        'reconstruction': 2,
        'tolerance': 1e-14,
        'ini': 'idefix.ini',
        'testfile': dir + '/test/pb2/testme.json',
        'testname': 'pb2'
      }, marks=(), id='pb2-idefix.ini'
    ),
    pytest.param(
      {
        'dumpname': 'dump.0001.dmp',
        'noplot': False,
        'reconstruction': 2,
        'tolerance': 1e-14,
        'ini': 'idefix-implicit.ini',
        'testfile': dir + '/test/pb2/testme.json',
        'testname': 'pb2'
      }, marks=(), id='pb2-idefix-implicit.ini'
    ),
  ]
