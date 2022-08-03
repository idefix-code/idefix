.. _contributing:

=======================
Contributing to Idefix
=======================

pre-commit framework
----------------------
Idefix is developed with the help of the `pre-commit <https://pre-commit.com>`_ framework.
We use `cpplint <https://en.wikipedia.org/wiki/Cpplint>`_ to validate code style, mostly
following the Google standards for C++, and several pre-commit hooks to automatically fix
some coding bad practices.

It is recommended (though not mandatory) to install pre-commit by running the following
script from the top level of the repo

.. code-block:: shell

  python3 -m pip install pre-commit
  pre-commit install


Then, as one checks in their contribution with ``git commit``, pre-commit hooks may perform
changes in situ. One then needs to re-add and enter the ``git commit`` command again for the
commit to be validated.
Note that an important hook that does *not* perform auto-fixes is ``cpplint``, so contributors
need to accomodate for this one by hand.

.. note::
  Note that if for any reason you do not wish, or are unable to install pre-commit in your
  environment, formatting errors will be caught by our CI after you open a merge-request.


Coding style guidelines
-----------------------

It is recommended to follow *Idefix* coding style standards whenever you develop your own setup
or physics module. These guidelines are mandatory if you intend to share your developements
in the *Idefix* git repository. Note that some of these guidelines will be enforced by the cpplint
included in ``pre-commit`` (see above).

#. Tabulation are two spaces
#. Class definition (.hpp) and implementation (.cpp) should be located in files with the same name as the class. For complex classes, one should have one directory
   with the class name, which contains at least the header (.hpp) and the class constructor implementation in files corresponding to the class names. One can then
   add the other class methods in other .cpp/.hpp files in the class directory, naming them after the methods they define. example:

   ::

     | MyComplexClass/
     |-- MyComplexClass.cpp
     |-- MyComplexClass.hpp
     |-- DoSomething.cpp
     |-- DoSomethingElse.hpp

#. Make your definitions and declaration as clear as possible (use words, not abbreviations!). Use camel cases
   if you want to use several word (e.g. ``mySmartVariable``).
#. Avoid underscore (_) and dash (-) as much as possible.
#. Class definitions should start with an upper case (e.g. ``DataBlock``).
#. Method names should start with an upper case. The first word of a method name is a verb (because methods do something). Example: ``MyClass::ComputeFibonacci``.
#. variable names and class instances should start with lower cases. (e.g. ``int myVariable``).
#. Macro and preprocessor functions should be all in upper case (e.g. ``EXPAND``)
#. If global variables//functions are needed, they should be in the namespace ``idfx::``. They follow the rules above for naming.


Running tests locally
---------------------

Tests are organized in categories such as hydro, MHD, MPI ... each corresponding to a bash script
located in ``$IDEFIX_DIR/test/``. For instance, the following invoke will run tests in batch and
without producing artifacts: compilation and output files are stored in temporary directories

.. code-block:: shell

   test/checks_hydro.sh


Requirements
============

Most tests are performed using Python scripts to analyse outputs.
It is recommended to use a Python virtual environment to run tests.
Tests scripts will look for an environment directory at ``$IDEFIX_DIR/.python-env/``.
To create this directory, run

.. code-block::shell

    cd $IDEFIX_DIR
    python -m venv .python-env
    source .python-env/bin/activate
    python -m pip install --upgrade pip
    python -m pip install -r test/python_requirements.txt
