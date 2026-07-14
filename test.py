#!/usr/bin/env python3

import os
import sys
import pytest

# set IDEFIX_DIR
source_dir = os.path.dirname(__file__)
os.environ["IDEFIX_DIR"] = source_dir

# idefix test class
sys.path.append(source_dir)
from pytools.idfx_test_run import IdexPytestRunner

# should be global so it remember state (last build) and avoids
# to build for each run if we just changed the ini file and run options.
gblIdefixPytestRunner = IdexPytestRunner(__file__)

# define the pytest test
@pytest.mark.parametrize("config", gblIdefixPytestRunner.genTests())
def test_idefix_build_run_check(config):
  gblIdefixPytestRunner.run(config)

# if called directly as a script
if __name__ == "__main__":
  gblIdefixPytestRunner.main(all=True)
