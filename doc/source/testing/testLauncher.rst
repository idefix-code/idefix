==========================
Test launcher and reporter
==========================

Overview
--------

The class :doc:`idfxTest <idfxTest>` provides the basement to implement an *Idefix* integration test for validation.
In order to ease launching all the tests, the user might prefer to use directly the ``./test.py`` command at the
root of the *Idefix* sources.

This script will run all the listed variants of *Idefix* and build a report in the terminal. At the
end of the run a standard ``junit.xml`` file is produced. This one can be translated into a browsable
HTML file.

Depencencies
------------

Before using :doc:`idfxTest <idfxTest>` you need to install some Python Depencencies (possibly in a ``virtual env``):

.. code-block:: shell

    pip install -r test/python_requirements.txt

Running
-------

To run the test you can basically :

.. code-block:: shell

    # Run all tests
    ./test.py

    # Run all tests in ./tests/HD
    ./test.py -subdir=./tests/HD

    # Select in more details the tests containing the "single" keyword
    # See pytest documentation for the exact advanced semantic
    ./test.py -subdir=./tests/HD -k single

    # Run in verbose
    ./test.py -v

The result of the execution will be an output like :

.. code-block:: text

    ============================================ test session starts =============================================
    collected 52 items / 44 deselected / 8 selected

    test.py::test_idefix_build_run_check[HD/sod-iso-idefix.ini-single-reconstruction-2] PASSED              [ 12%]
    test.py::test_idefix_build_run_check[HD/sod-iso-idefix-hll.ini-single-reconstruction-2] PASSED          [ 25%]
    test.py::test_idefix_build_run_check[HD/sod-iso-idefix-hllc.ini-single-reconstruction-2] PASSED         [ 37%]
    test.py::test_idefix_build_run_check[HD/sod-iso-idefix-tvdlf.ini-single-reconstruction-2] PASSED        [ 50%]
    test.py::test_idefix_build_run_check[HD/sod-idefix.ini-single-reconstruction-2] PASSED                  [ 62%]
    test.py::test_idefix_build_run_check[HD/sod-idefix-hll.ini-single-reconstruction-2] PASSED              [ 75%]
    test.py::test_idefix_build_run_check[HD/sod-idefix-hllc.ini-single-reconstruction-2] PASSED             [ 87%]
    test.py::test_idefix_build_run_check[HD/sod-idefix-tvdlf.ini-single-reconstruction-2] PASSED            [100%]

    ---------- generated xml file: idefix-tests.junit.xml ---------------------------------------------------------
    =============================== 8 passed, 44 deselected in 73.03s (0:01:13) ===================================

When seeing an error, the output of the command will be printed at the end.

HTML report
-----------

If you want to to generate an HTML page from the report you can proceed by using the Python package ``junit2html`` :

.. code-block:: shell

    # install junit2html
    pip install junit2html

    # convert the report
    junit2html ./idefix-tests.junit.xml ./idefix-tests.junit.html

Advanced usage of the command
-----------------------------

Here the options supported by the test script :

.. code-block:: text

    usage: test.py [-h] [-noplot] [-ploterr] [-cmake CMAKE [CMAKE ...]] [-definitions DEFINITIONS]
                [-dec DEC [DEC ...]] [-check] [-cuda] [-intel] [-hip] [-single] [-vectPot]
                [-reconstruction RECONSTRUCTION] [-idefixDir IDEFIXDIR] [-mpi] [-all] [-init]
                [-Werror] [-ccache] [-restart] [-v] [--help-pytest] [-fake] [-subdir SUBDIR]

    options:
    -h, --help            show this help message and exit
    -noplot               disable plotting in standard tests
    -ploterr              Enable plotting on error in regression tests
    -cmake CMAKE [CMAKE ...]
                            CMake options
    -definitions DEFINITIONS
                            definitions.hpp file
    -dec DEC [DEC ...]    MPI domain decomposition
    -check                Only perform regression tests without compilation
    -cuda                 Test on Nvidia GPU using CUDA
    -intel                Test compiling with Intel OneAPI
    -hip                  Test on AMD GPU using HIP
    -single               Enable single precision
    -vectPot              Enable vector potential formulation
    -reconstruction RECONSTRUCTION
                            set reconstruction scheme (2=PLM, 3=LimO3, 4=PPM)
    -idefixDir IDEFIXDIR  Set directory for idefix source files (default $IDEFIX_DIR)
    -mpi                  Enable MPI
    -all                  Do all test suite (otherwise, just do the test with the current configuration)
    -init                 Reinit reference files for non-regression tests (dangerous!)
    -Werror               Consider warnings as errors
    -ccache               Use ccache to reduce the build time over multiple run of the test suite.
    -restart              Enable creating a restart from a checkpoint.
    -v, --verbose         Enable verbose mode, by not capturing the output.
    --help-pytest         Display the options you can transmit directly to pytest in addition to the specific to idefix tests.
    -fake                 Make a fake run by just logging the actions to validate that we generate same command over refactoring.
    -subdir SUBDIR        Select the test in the given subdir not to run all.

The script is built on top of the ``pytest`` command so you automatically get access
to all the advanced option this command provide. Here a few examples :

.. code-block:: shell

    # get the pytest help
    ./test.py --help-pytest

    # let pytest filtering the tests
    ./test.py -k "single and mpi"

    # stop on first failure
    ./test.py -x

    # re-run only the last failed tests
    ./test.py --lf

Definition of the tests
-----------------------

The script is simply searching all the files names ``testme.json`` into the ``test`` directory.
This file describe the combination of parameters to use to produce the list of *Idefix* run to perform.
For a single basic configuration one can use :

.. code-block:: json

    {
        "variants": {
            "dumpname": "dump.0001.dmp",
            "ini": "idefix.ini",
            "vectPot": false,
            "single": false,
            "reconstruction": 2,
            "mpi": false,
            "standardTest": false,
            "tolerance": 0
        }
    }

Available parameters
--------------------

The parameters in the ``variants`` dictionnary correspond to the options supported by the
:doc:`idfxTest <idfxTest>` script to configure the build and run of *Idefix*.

In addition there is some extra keys which are dedicated to the json interpretation layer :

.. list-table::
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``dumname``
     - ``dump.0001.dmp``
     - Dump file to use to check the results after the run.
   * - ``ini``
     - ``idefix.ini``
     - The configuration file to use.
   * - ``tolerance``
     - ``0``
     - The margins to allow when checking the results.
   * - ``standardTest``
     - ``true``
     - Runs any Python-based standard tests (e.g., ``testidefix.py``) present in the test directory for additional validation.
   * - ``nonRegressionTest``
     - ``true``
     - Compares the output dump file to a reference file using RMSE; fails if the error exceeds the tolerance.
   * - ``nonRegressionTestIni``
     - Same than ``ini``
     - When making restart you might want to make the check using the inirial configuration file.
   * - ``multirun``
     - ``{}``
     - See the multi-run section below.

Looping over parameters
-----------------------

You might want to explore running Idefix within parameter ranges (configuration files, modes).
For this simply list the values you want as a list. The test script will automatically
generate all combinations.

.. code-block:: json

    {
        "variants": {
            "dumpname": "dump.0001.dmp",
            "ini": ["idefix.ini","idefix-hll.ini"],
            "vectPot": [false, true],
            "single": false,
            "reconstruction": 2,
            "mpi": [false, true],
            "standardTest": false,
            "tolerance": 0
        }
    }

It will automatically produce the tests :

* HD/sod-iso-idefix.ini
* HD/sod-iso-idefix.ini-vectPot
* HD/sod-iso-idefix.ini-mpi
* HD/sod-iso-idefix.ini-vectPot-mpi
* HD/sod-iso-idefix-hll.ini
* HD/sod-iso-idefix-hll.ini-vectPot
* HD/sod-iso-idefix-hll.ini-mpi
* HD/sod-iso-idefix-hll.ini-vectPot-mpi

Specific keys
-------------

There is some keys which are by default some arrays, they will not be considered
as combination rules :

* ``dec``
* ``multirun``
* ``restart_no_overwrite``
* ``tolerance``

Reduce the combinations
-----------------------

You might not want to see all the combinations but just a few, for this, you
can list several sets as a list. Here using single only on half of the modes.

.. code-block:: json

    {
        "variants": [
            {
                "dumpname": "dump.0001.dmp",
                "ini": ["idefix.ini"],
                "vectPot": false,
                "single": false,
                "reconstruction": 2,
                "mpi": [false, true],
                "standardTest": false,
                "tolerance": 0
            },{
                "dumpname": "dump.0001.dmp",
                "ini": ["idefix-hll.ini"],
                "vectPot": true,
                "single": true,
                "reconstruction": 2,
                "mpi": [false, true],
                "standardTest": false,
                "tolerance": 0
            }
        ]
    }

* HD/sod-iso-idefix.ini
* HD/sod-iso-idefix.ini-mpi
* HD/sod-iso-idefix-hll.ini-single-vectPot
* HD/sod-iso-idefix-hll.ini-single-vectPot-mpi

Naming the test
---------------

If you prefer to see the options appearing in a specific order in the generated test name,
you can provide the key ``namings`` listing as a comma separated list the variables in
the order you want to see them composing the test name.

.. code-block:: json

    {
        "namings": "ini,single,mpi",
        "variants": {
            "dumpname": "dump.0001.dmp",
            "ini": ["idefix.ini","idefix-hll.ini"],
            "vectPot": false,
            "single": [false, true],
            "reconstruction": 2,
            "mpi": [false, true],
            "standardTest": false,
            "tolerance": 0
        }
    }


By default, the alphabetical order will be used.

When clauses
------------

You can also dynamically override a specific value when a parameter value is selected.
It is just like if you used and IF statement.

.. code-block:: json

    {
        "namings": "ini,single,mpi",
        "variants": {
            "dumpname": "dump.0001.dmp",
            "ini": ["idefix.ini","idefix-hll.ini"],
            "vectPot": false,
            "single": [false, true],
            "reconstruction": 2,
            "mpi": [false, true],
            "standardTest": false,
            "tolerance": 0
        },
        "when": {
            "conditions": {
                "single": true
            },
            "apply": {
                "reconstruction": 1
            }
        }
    }

In this case, ``reconstruction`` will be set to 1 when ``single`` is equal to ``true``.

Note that the ``conditions`` field can contains several values which will be treaded
as and AND logical operator.

You can provide several ``when`` clauses by using a list of them instead of directly
the dictionnary.

Making multi-run steps
----------------------

In order to validate checkpoint restart, or continuing a simulation with dirfferent
tunnings the script supports decribing multi-run configurations.

They are described like :

.. code-block:: json

    {
        "variants": {
            "dumpname": "dump.0001.dmp",
            "ini": ["idefix.ini"],
            "vectPot": false,
            "single": false,
            "reconstruction": 2,
            "mpi": false,
            "standardTest": false,
            "tolerance": 0,
            "dec": [2,2,2],
            "multirun": [
                {
                },{
                    "mpi": true,
                    "restart": true,
                    "restart_no_overwrite": ["dump.0001.dmp", "data.0005.vtk"]
                }
            ]
        },
    }

Using the idfxTest options
--------------------------

As the ``test.py`` uses the class described in :doc:`idfxTest <idfxTest>` it also supports all the
command line options it offers.

Most usefull might be the enabling of ``ccache`` to reduce the compilation time from one
run to another.

.. code-block:: shell

    ./test.py -ccache
