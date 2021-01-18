=====================
Reference guide
=====================
First, it is good practice to set the environment variable ``IDEFIX_DIR`` to the root path of your *Idefix* distribution, as it is needed at several stages. Setting up *Idefix* for a particular problem implies editing several files, some of which are automatically generated. There are essentially 4 files to a proper *Idefix* setup:

``definitions.hpp``
    The problem header file. It contains C preprocessor directives for configuration options which require a re-compilation of the code.

``setup.cpp``
    The physical setup definition. Its role is to define the initial conditions of the problem, but also to provide functions for user-defined
    properties (boundary conditions, diffusivity, potential, etc.)

``idefix.ini``
    The input file. This file is read at runtime *only*. It defines most of the properties of the code: resolution, integrator, output, physical modules.

``Makefile``
    The compilation script. This file is automatically created during the code configuration by the configure script. Although it is possible to edit and modify
    the makefile, it is recommended to change options when calling the configure script instead.

Usually, one starts a new project by copying an example from the test suite of idefix in ``$IDEFIX_DIR/test``. From there, one modifies the setup for the particular problem at hand.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   reference/definitions.hpp
   reference/makefile
   reference/setup.cpp
   reference/idefix.ini
   reference/commandline









Migrating from PLUTO
====================
