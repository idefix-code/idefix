#####################################################################################
# Idefix MHD astrophysical code
# Copyright(C) Sébastien Valat <sebastien.valat@univ-grenoble-alpes.fr>
# and other code contributors
# Licensed under CeCILL 2.1 License, see COPYING for more information
#####################################################################################

import copy
import pytest

DO_NOT_LOOP_ON = ['restart_no_overwrite', "dec", "multirun", "check_file_produced"]

class IdefixDirTestGenerator:
  '''
  Class used to generate the various configuration to run by parsing the
  files `testme.py` found in the hierarchy of the /test directory of Idefix.
  '''

  def __init__(self, currentTestFile: str, name: str = ""):
    '''
    Constructor or the class.

    Args:
      currentTestFile (str):
        Path the the current python file. Normally you
        simply pass __file__ to this parameter.
      name (str):
        Define a name for the test, can be empty.
    '''

    self.currentTestFile = currentTestFile
    self.currentTestName = name

  # generate the list of configs to run
  def genTestConfigs(self, names:str, params, whenClauses = {}) -> list:
    '''
    Generate the the list of configurations as pytest parameters.
    It will unpack the configuration set by looping on all combinations defined
    by the given sets.

    Args:
      names (str):
        Comma separated list of variables do consider to build the
        name of the file.
      params (dict|list):
        A configuration set as a dictionnary or a list of
        configuration set.
      whenCaluses (dict|list):
        Provide a set of clauses to apply after unpacking
        the configuration so we can patch some values depending on some others.

    Returns:
      A list of pytest.param() ready to be fiven to parametrized pytest functions.
    '''
    # get name ordering list
    nameList = names.split(',')

    # gen list of complete configs
    all_configs = []
    if isinstance(params, dict):
      all_configs += self._genOneConfigSeries(names, params)
    elif isinstance(params, list):
      for p in params:
        all_configs += self._genOneConfigSeries(names, p)
    else:
      raise Exception("Should never be called !")

    # convert as parametrize with nice name
    result = []
    for config in all_configs:
      # append the file
      config['testfile'] = self.currentTestFile
      config['testname'] = self.currentTestName
      # gen name
      nameParts = [self.currentTestName]
      for name in nameList:
        if isinstance(config[name], bool):
          if config[name]:
            nameParts.append(name)
        elif isinstance(config[name], str):
          nameParts.append(str(config[name]))
        else:
          nameParts.append(f"{name}-{config[name]}")
      confName = "-".join(nameParts)

      # apply when clause
      config = self._applyWhen(config, whenClauses)

      result.append(pytest.param(config, id=confName))

    # ok
    return result

  def extractNamingParameters(self, params) -> list:
    '''
    Loop on the parameters and check automatically what are the list of variable parameters.

    Args:
      params (list|dict):
        The list of configuration sets to scan or a single set.

    Returns:
      The list of names of the variable parameters.
    '''

    # if not a list make a list
    if not isinstance(params, list):
      params = [params]

    # see params
    seen = {}
    variables = []

    # loop
    for param_set in params:
      for key, value in param_set.items():
        if key in variables:
          pass
        elif key in DO_NOT_LOOP_ON:
          pass
        elif isinstance(value, list):
          variables.append(key)
        elif key in seen and seen[key] != value:
          variables.append(key)
        elif key not in seen:
          seen[key] = value

    # by default sort by alphabetic order about var names
    # TODO make something better by assigning a priority to vars
    variables.sort()

    # ok
    return ','.join(variables)

  def _genNextLevelCombinations(self, input: list, paramName: str, paramValues: list) -> list:
    '''
    Take an input case list and unpack the given parameter to build the new combinations.

    Args:
      input (list):
        The incoming list of combinations already unpacked before this call.
      paramName (str):
        Name of the parameter to unpack.
      paramValues (list):
        The list of values to unpack and to build combinations for.

    Returns:
      The updated list of run sets.
    '''
    result = []
    for entry in input:
      for value in paramValues:
        v = copy.deepcopy(entry)
        v[paramName] = value
        result.append(v)
    return result

  def _genOneConfigSeries(self, names: str, config: dict) -> list:
    '''
    Generate the the list of configurations as pytest parameters.
    It will unpack the configuration set by looping on all combinations defined
    by the given sets.

    Args:
      names (str):
        Comma separated list of variables do consider to build the
        name of the file.
      config (dict):
        A configuration set as a dictionnary.

    Returns:
      A list of pytest.param() ready to be fiven to parametrized pytest functions.
    '''
    # get name ordering list
    nameList = names.split(',')

    # if there is ini in the list we put it at the end
    loopOrder = nameList.copy()
    if 'ini' in loopOrder:
      loopOrder.remove('ini')
      loopOrder.append('ini')

    # init core with everything not a list
    core = {}
    for key, value in config.items():
      if isinstance(value, list) and key not in DO_NOT_LOOP_ON:
        assert key in nameList, f"All variable parameteres should be ordered in the names list, '{key}' is not."
      else:
        core[key] = copy.deepcopy(value)

    # at start we have only default core
    result = [core]

    # loop
    for key in loopOrder:
      value = config[key]
      assert isinstance(value, list), f"This parameter is marked as a list but is not a list : {key}={value} !"
      result = self._genNextLevelCombinations(result, key, value)

    # ok
    return result

  def _matchWhenClause(self, config: dict, when_clause: dict) -> bool:
    '''
    Check the matching of a given when clause on the given configuration.

    Args:
      config (dict): The configuration.
      when_clause (dict): The when clause to check.

    Returns:
      True if the clause, match, False otherwise.
    '''
    for key, value in when_clause.items():
      if config[key] != value:
        return False
    return True

  def _applyWhen(self, config: dict, when: dict) -> dict:
    '''
    Check if the when clause applies and apply if if any.

    Args:
      config (dict): The configuration.
      when_clause (dict): The when clause to check and apply.

    Returns:
      The fixed config.
    '''
    # nothing to do
    if when == {}:
      return config

    # clone
    result = copy.deepcopy(config)

    # loop on when
    if isinstance(when, list):
      for clause in when:
        if self._matchWhenClause(config, clause['conditions']):
          result.update(clause['apply'])
    elif isinstance(when, dict):
      if self._matchWhenClause(config, when['conditions']):
        result.update(when['apply'])
    else:
      raise Exception(f"Invalid type for 'when' : {when}")

    # ok
    return result
