=========
idfxTest
=========

.. autoclass:: idfxTest
   :members:
   :undoc-members:
   :show-inheritance:

Overview
--------

The ``idfxTest`` class provides a high-level interface for automating the configuration, compilation, execution, and regression testing of Idefix simulations. It is designed to be used in test scripts (such as ``testme.py``) to streamline the testing workflow, including handling reference files and plotting differences.

Constructor and Command-Line Options
------------------------------------

The constructor parses command-line arguments using ``argparse``. These options can be passed directly to the test script or via the command line. The following options are available:

.. list-table::
   :header-rows: 1

   * - Option
     - Attribute
     - Description
   * - ``-noplot``
     - ``noplot``
     - Disable plotting in standard tests (default: True).
   * - ``-ploterr``
     - ``ploterr``
     - Enable plotting of differences when regression tests fail.
   * - ``-cmake OPT [OPT ...]``
     - ``cmake``
     - Extra CMake options (list of strings).
   * - ``-definitions FILE``
     - ``definitions``
     - Specify a custom ``definitions.hpp`` file.
   * - ``-dec NX NY NZ``
     - ``dec``
     - MPI domain decomposition (list of integers).
   * - ``-check``
     - ``check``
     - Only perform regression tests without compilation.
   * - ``-cuda``
     - ``cuda``
     - Enable CUDA backend for Nvidia GPUs.
   * - ``-intel``
     - ``intel``
     - Use Intel OneAPI compilers.
   * - ``-hip``
     - ``hip``
     - Enable HIP backend for AMD GPUs.
   * - ``-single``
     - ``single``
     - Enable single precision.
   * - ``-vectPot``
     - ``vectPot``
     - Enable vector potential formulation.
   * - ``-reconstruction N``
     - ``reconstruction``
     - Set reconstruction scheme (2=PLM, 3=LimO3, 4=PPM).
   * - ``-idefixDir PATH``
     - ``idefixDir``
     - Set the directory for Idefix source files (default: ``$IDEFIX_DIR``).
   * - ``-mpi``
     - ``mpi``
     - Enable MPI parallelism.
   * - ``-all``
     - ``all``
     - Run the full test suite with multiple configurations.
   * - ``-init``
     - ``init``
     - Reinitialize reference files for non-regression tests.
   * - ``-Werror``
     - ``Werror``
     - Treat compiler warnings as errors.

Main Methods
------------

.. list-table::
   :header-rows: 1

   * - Method
     - Description
   * - ``configure``
     - Runs CMake to configure the build system for Idefix, using options set by command-line flags (e.g., precision, MPI, CUDA, etc.).
   * - ``compile``
     - Compiles the Idefix code using ``make`` with the specified number of parallel jobs.
   * - ``run``
     - Executes the Idefix binary, optionally with MPI, using the provided input file and runtime options.
   * - ``checkOnly``
     - Performs regression testing only, without compiling or running the code (useful for checking outputs after a manual run).
   * - ``standardTest``
     - Runs any Python-based standard tests (e.g., ``testidefix.py``) present in the test directory for additional validation.
   * - ``nonRegressionTest``
     - Compares the output dump file to a reference file using RMSE; fails if the error exceeds the tolerance.
   * - ``compareDump``
     - Compares two arbitrary dump files using the same logic as ``nonRegressionTest``.
   * - ``makeReference``
     - Copies the specified output file to the reference directory, updating the reference for future regression tests.

Usage Example
-------------

Below is an example inspired by ``testme.py`` from ``test/HD/sod/testme.py``. This demonstrates a typical workflow for running tests and performing regression checks.

.. code-block:: python

   import pytools.idfx_test as tst

   name = "dump.0001.dmp"

   def testMe(test):
       test.configure()
       test.compile()
       inifiles = ["idefix.ini", "idefix-hll.ini", "idefix-hllc.ini", "idefix-tvdlf.ini"]
       if test.reconstruction == 4:
           inifiles = ["idefix-rk3.ini", "idefix-hllc-rk3.ini"]

       # Loop over all ini files for this test
       for ini in inifiles:
           test.run(inputFile=ini)
           if test.init:
               test.makeReference(filename=name)
           test.standardTest()
           test.nonRegressionTest(filename=name)

   test = tst.idfxTest()

   if not test.all:
       if test.check:
           test.checkOnly(filename=name)
       else:
           testMe(test)
   else:
       test.noplot = True
       for rec in range(2, 5):
           test.vectPot = False
           test.single = False
           test.reconstruction = rec
           test.mpi = False
           testMe(test)

       # Test in single precision
       test.reconstruction = 2
       test.single = True
       testMe(test)

How to Run
----------

You can run the test script from the command line, passing any of the supported options. For example:

.. code-block:: bash

   python testme.py -mpi -dec 2 2 -reconstruction 3 -single -ploterr -idefixDir /path/to/idefix

This will configure, compile, and run the test in MPI mode with a 2x2 domain decomposition, third-order reconstruction, single precision, and plotting enabled for regression errors.

Reference File Management
------------------------

- Reference files are stored in ``$IDEFIX_DIR/reference/<testDir>``.
- The filename is generated based on precision, reconstruction, input file, and vector potential settings.
- Use ``test.init`` to regenerate reference files (dangerous: overwrites existing references).

Regression Testing
------------------

- The ``nonRegressionTest`` method compares output dumps to reference files using RMSE.
- If the error exceeds the tolerance, the test fails and (optionally) plots the difference.
