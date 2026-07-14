#####################################################################################
# Idefix MHD astrophysical code
# Copyright(C) Sébastien Valat <sebastien.valat@univ-grenoble-alpes.fr>
# and other code contributors
# Licensed under CeCILL 2.1 License, see COPYING for more information
#####################################################################################

import os
import sys
import json
import glob
import copy
import pytest
# idefix test class
import pytools.idfx_test as tst
from  pytools.idfx_test_gen import IdefixDirTestGenerator
from contextlib import contextmanager

@contextmanager
def moveInDir(path):
    '''
    Change current working dir for the given directoy for what is inside the
    with statement.
    '''
    oldpwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(oldpwd)

class IdexPytestRunner:
  '''
  Implement the Idefix pytest runner to scan and run all the tests described by the files
  testme.json into the /test directory of Idefix sources.
  '''

  def __init__(self, parentScriptFile: str):
    self.currentTestFile=""
    self.currentTestRunner: tst.idfxTest=None
    self.parentScritFile=parentScriptFile
    self.filterSubdir = os.environ.get("IDEFIX_TEST_FILTER_SUBDIR", "./test/")

  def _validateNaming(self, namings: str, autoExtracted: str, file: str):
    names = namings.split(',')
    for n in autoExtracted.split(','):
      if not n in names:
        raise Exception(f"Naming parameter list not match the auto-detectection, some are missinge. You gave : '{names}', detected '{autoExtracted}' into {file}")

  def _makeVariableArgAsList(self, namings: str, variants) -> list:
    # make as a list
    if isinstance(variants, list) == False:
      variants = [variants]

    # split namings
    namings = namings.split(",")

    # loop on each
    for variant in variants:
      for name in namings:
        if name in variant and isinstance(variant[name], list) == False:
          variant[name] = [variant[name]]

  def genTests(self) -> list:
    sourceDir = os.path.dirname(self.parentScritFile)

    # loop over all tests
    result = []
    with moveInDir(sourceDir):
      # if missing /
      if self.filterSubdir != "" and self.filterSubdir[-1] != "/":
        self.filterSubdir += "/"

      # walk in test dir to find the tests & sort by name
      testfiles = glob.glob(self.filterSubdir + "**/testme.json", recursive=True)
      testfiles.sort()

      # loop over each
      for testfile in testfiles:
        try:
          # calc some paths
          testfileRelPath = os.path.relpath(testfile, os.path.join(sourceDir, 'test'))
          testfilePath = os.path.abspath(os.path.join('test',testfileRelPath))
          testfileDir = os.path.dirname(testfileRelPath)

          # load json & build the inner test combinations
          with open(testfilePath, 'r') as fp:
            test = json.load(fp)
            idefixTestGenerator=IdefixDirTestGenerator(testfilePath, testfileDir)
            if 'namings' in test:
              namings = test['namings']
              autoExtracted = idefixTestGenerator.extractNamingParameters(test['variants'])
              self._validateNaming(namings, autoExtracted, testfilePath)
            else:
              namings = idefixTestGenerator.extractNamingParameters(test['variants'])

            # required to simplify the algos later, if var is listed as variable, we need to loop over it.
            self._makeVariableArgAsList(namings, test['variants'])

            # gen
            result += idefixTestGenerator.genTestConfigs(namings, test['variants'], test.get('when', {}))
        except Exception as e:
          raise Exception(f"Fail to generate tests from {testfileRelPath} : {e}")

    # ok
    return result

  def run(self, config: dict) -> None:
    # clone before modify to not modity for caller
    config = copy.deepcopy(config)

    # print config
    print("***************************************************")
    print(json.dumps(config, indent='\t'))
    print("***************************************************")

    # extract some infos for local usage
    testfile = config["testfile"]
    dumpname = config['dumpname']
    testname = config["testname"]
    tolerance = config.get("tolerance", 0)
    definitionFile = config.get("definitionFile", "")
    standardTest = config.get("standardTest", True)
    nonRegressionTest = config.get("nonRegressionTest", True)
    nonRegressionTestIni = config.get("nonRegressionTestIni", None)
    check_file_produced = config.get("check_file_produced", [])
    problemDir = os.path.dirname(testfile)

    # cleanup some keyword not handled at the
    # level of idx_test so we don't perturbate it
    del config['dumpname']
    if 'definitionFile' in config:
      del config['definitionFile']
    if 'tolerance' in config:
      del config['tolerance']
    if 'standardTest' in config:
      del config['standardTest']
    if 'nonRegressionTest' in config:
      del config['nonRegressionTest']
    if 'nonRegressionTestIni' in config:
      del config['nonRegressionTestIni']

    # if switch from test, rebuild the runner (a runner make for one dir)
    if self.currentTestFile != testfile:
      self.currentTestRunner = tst.idfxTest(testfile, name=testname)
      self.currentTestFile = testfile

    # run
    with moveInDir(problemDir):
      self._runNonRegression(dumpname, config['ini'], config, tolerance=tolerance, definitionFile=definitionFile, standardTest=standardTest, nonReg=nonRegressionTest, nonRegIni=nonRegressionTestIni)

    # check produced
    for file in check_file_produced:
      if not os.path.exists(file):
        raise Exception(f"Don't find expected file to be produced by the run : {file} !")

  def _runNonRegression(self, dumpname, ini, config_override, tolerance=0, definitionFile="", nonReg=True, nonRegIni=None, standardTest=True, first_run_ini=None,first_run_dumpname=None,configure_and_compile=True):
    if 'multirun' in config_override:
      self._runNonRegMultirun(dumpname, ini, config_override, tolerance=tolerance, nonReg=nonReg, nonRegIni=nonRegIni, standardTest=standardTest, configure_and_compile=configure_and_compile, definitionFile=definitionFile, first_run_ini=first_run_ini, first_run_dumpname=first_run_dumpname)
    else:
      # single basic run
      self._runNonRegSingleRun(dumpname, ini, config_override, tolerance=tolerance, nonReg=nonReg, standardTest=standardTest, configure_and_compile=configure_and_compile, definitionFile=definitionFile, first_run_ini=first_run_ini, first_run_dumpname=first_run_dumpname)

  def _runNonRegMultirun(self, dumpname, ini, config_override, tolerance=0, definitionFile="", nonReg=True, nonRegIni=None, standardTest=True, first_run_ini=None,first_run_dumpname=None,configure_and_compile=True):
    # check
    assert 'multirun' in config_override

    # loop over runs
    for run in config_override['multirun']:
      # copy config
      run_config = copy.deepcopy(config_override)

      # patch a bit
      del run_config['multirun']
      run_config.update(run)
      nonReg=run_config.get('nonRegressionTest', nonReg)
      dumpname=run_config.get('dumpname', dumpname)
      if 'nonRegressionTest' in run_config:
        del run_config['nonRegressionTest']
      standardTest=run_config.get('standardTest', standardTest)
      if 'standardTest' in run_config:
        del run_config['standardTest']

      # make single run
      self._runNonRegSingleRun(dumpname, run_config['ini'], run_config, definitionFile=definitionFile, tolerance=tolerance, nonReg=nonReg, nonRegIni=nonRegIni, standardTest=standardTest)

  def _runNonRegSingleRun(self, dumpname, ini, config_override, tolerance=0, definitionFile="", nonReg=True, nonRegIni=None, standardTest=True, first_run_ini=None,first_run_dumpname=None,configure_and_compile=True):
    # build the runner
    idefixTest = self.currentTestRunner

    # handle special override which should not go into
    # idefixTest because it is not supported by the idfx_test layer.
    config_override = copy.deepcopy(config_override)
    nonRegIni = config_override.get("nonRegressionTestIni", nonRegIni)
    if 'nonRegressionTestIni' in config_override:
      del config_override['nonRegressionTestIni']

    # apply config
    idefixTest.applyConfig(config_override)

    # recompile if needed
    if configure_and_compile:
      idefixTest.configure(override=config_override, definitionFile=definitionFile)
      idefixTest.compile()

    if first_run_ini:
      idefixTest.run(inputFile=first_run_ini)
      if nonReg:
        idefixTest.nonRegressionTest(filename=first_run_dumpname, tolerance=tolerance)

    # Test the restart option
    file_mtime={}
    if not idefixTest.fake:
      for file in idefixTest.restart_no_overwrite:
        file_mtime[file] = os.path.getmtime(file)

    # restart
    if idefixTest.restart:
      restart = 1
    else:
      restart = -1

    # run
    idefixTest.run(inputFile=ini, restart=restart)

    # regen ref if needed
    if idefixTest.init:
      idefixTest.makeReference(filename=dumpname)

    # check outputs
    if standardTest:
      idefixTest.standardTest()
    if nonReg:
      if nonRegIni:
        idefixTest.inifile = nonRegIni
      idefixTest.nonRegressionTest(filename=dumpname, tolerance=tolerance)

    # check that we didn't overrite the file during the restart
    if not idefixTest.fake:
      for file in idefixTest.restart_no_overwrite:
        assert file_mtime[file] == os.path.getmtime(file), f"Dump file {file} was overwritten on restart"

  def main(self, all: bool = False):
    if all:
      sys.argv.append("-all")
    idefixTest = tst.idfxTest(self.parentScritFile, name="main")
    os.environ["IDEFIX_TEST_FILTER_SUBDIR"] = idefixTest.filterSubdir

    if idefixTest.all:
      pytest.main(['-v', '--no-header', '--junit-xml=idefix-tests.junit.xml', '--tb=short'] + idefixTest.remainingArgs + [self.parentScritFile])
    else:
      assert False, "Not yet supported !"
    #elif self.check:
    #  idefixTest.checkOnly(filename=dumpname, tolerance=tolerance)
    #else:
    #  for ini in ini_list:
    #      self.runNonRegression(dumpname, ini, {}, tolerance=tolerance)
