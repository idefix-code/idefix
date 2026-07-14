#####################################################################################
# Idefix MHD astrophysical code
# Copyright(C) Sébastien Valat <sebastien.valat@univ-grenoble-alpes.fr>
# and other code contributors
# Licensed under CeCILL 2.1 License, see COPYING for more information
#####################################################################################

from ..idfx_test_gen import IdefixDirTestGenerator
import pytest

def test_extractNamingParameters():
  # build generator
  gen = IdefixDirTestGenerator(__file__, "unit-test")

  lst = gen.extractNamingParameters([
    {
      "dec": [1,2,3],
      "a": 10,
      "b": 11,
      "c": [True, False]
    },
    {
      "dec": [1,2,4],
      "a": 10,
      "b": 12,
      "d": [True, False]
    },
  ])

  # valid
  assert lst == 'b,c,d'

def test_genNextLevelCombinations():
  # build generator
  gen = IdefixDirTestGenerator(__file__, "unit-test")

  # init with single element
  core = {}
  init = [core]

  # gen first set of combinations
  lst1 = gen._genNextLevelCombinations(init, "mpi", [True, False])
  assert lst1 == [
    {"mpi": True},
    {"mpi": False},
  ]

  # gen first second set of combinations
  lst3 = gen._genNextLevelCombinations(lst1, "name", ["a", "b"])
  assert lst3 == [
    {"mpi": True, "name": "a"},
    {"mpi": True, "name": "b"},
    {"mpi": False, "name": "a"},
    {"mpi": False, "name": "b"},
  ]

def test_genOneConfigSeries():
  # build generator
  gen = IdefixDirTestGenerator(__file__, "unit-test")

  # gen
  result = gen._genOneConfigSeries("mpi,name", {
    "mpi": [True, False],
    "name": ["a", "b"]
  })

  # check
  assert result == [
    {
      'mpi': True,
      'name': 'a',
    },{
      'mpi': True,
      'name': 'b',
    },{
      'mpi': False,
      'name': 'a',
    },{
      'mpi': False,
      'name': 'b',
    },
  ]

def test_matchWhenClause_single():
  # build generator
  gen = IdefixDirTestGenerator(__file__, "unit-test")

  # calls
  assert gen._matchWhenClause({"a":10, "b": 11}, {"a":10}) == True
  assert gen._matchWhenClause({"a":10, "b": 11}, {"a":11}) == False
  assert gen._matchWhenClause({"a":10, "b": 11}, {"b":11}) == True
  assert gen._matchWhenClause({"a":10, "b": 11}, {"b":10}) == False

def test_matchWhenClause_and():
  # build generator
  gen = IdefixDirTestGenerator(__file__, "unit-test")

  # calls
  assert gen._matchWhenClause({"a":10, "b": 11}, {"a":10, "b":11}) == True
  assert gen._matchWhenClause({"a":10, "b": 11}, {"a":11, "b":11}) == False

def test_applyWhen_none():
  # build generator
  gen = IdefixDirTestGenerator(__file__, "unit-test")

  # call when clause
  res = gen._applyWhen({
    "a": 10,
    "b": 11
  }, when={})

  # check result
  assert res == {"a": 10, "b": 11}

def test_applyWhen_single_apply_yes():
  # build generator
  gen = IdefixDirTestGenerator(__file__, "unit-test")

  # call when clause
  res = gen._applyWhen({
    "a": 10,
    "b": 11
  }, when={
    "conditions": {
      "a": 10,
    },
    "apply": {
      "b": 12
    }
  })

  # check result
  assert res == {"a": 10, "b": 12}

def test_applyWhen_single_apply_no():
  # build generator
  gen = IdefixDirTestGenerator(__file__, "unit-test")

  # call when clause
  res = gen._applyWhen({
    "a": 10,
    "b": 11
  }, when={
    "conditions": {
      "a": 9,
    },
    "apply": {
      "b": 12
    }
  })

  # check result
  assert res == {"a": 10, "b": 11}

def test_applyWhen_list_apply_yes():
  # build generator
  gen = IdefixDirTestGenerator(__file__, "unit-test")

  # call when clause
  res = gen._applyWhen({
    "a": 10,
    "b": 11
  }, when=[
    {
      "conditions": {
        "a": 10,
      },
      "apply": {
        "b": 12
      }
    },
    {
      "conditions": {
        "a": 13,
      },
      "apply": {
        "b": 13
      }
    }
  ])

  # check result
  assert res == {"a": 10, "b": 12}

def test_gen_full():
  # build generator
  gen = IdefixDirTestGenerator(__file__, "unit-test")

  # gen
  result = gen.genTestConfigs("mpi,name", {
    "mpi": [True, False],
    "name": ["a", "b"]
  }, {
    "conditions": {
      "mpi": True
    },
    "apply": {
      "extra": 10,
    }
  })

  # check
  assert result == [
    pytest.param({'mpi': True, 'name': 'a', 'testfile': __file__, 'testname': 'unit-test', 'extra': 10}, id='unit-test-mpi-a'),
    pytest.param({'mpi': True, 'name': 'b', 'testfile': __file__, 'testname': 'unit-test', 'extra': 10}, id='unit-test-mpi-b'),
    pytest.param({'mpi': False, 'name': 'a', 'testfile': __file__, 'testname': 'unit-test'}, id='unit-test-a'),
    pytest.param({'mpi': False, 'name': 'b', 'testfile': __file__, 'testname': 'unit-test'}, id='unit-test-b'),
  ]
