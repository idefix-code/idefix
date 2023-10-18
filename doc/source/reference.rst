=====================
User guide
=====================
First, it is good practice to set the environment variable ``IDEFIX_DIR`` to the root path of your *Idefix* distribution, as it is needed at several stages. Setting up *Idefix* for a particular problem implies editing several files, some of which are automatically generated.
To set up a particular *Idefix* setup, 5 steps are required.

1. Create/edit ``definitions.hpp``, the problem header file. It contains C preprocessor directives to set the code dimensionality, geometry and equation of state.

2. Create/edit ``setup.cpp``, the physical setup of your particular problem. Its role is to define the initial conditions of the problem, but also to provide user-defined functions (boundary conditions, diffusivity, potential, etc.) which will be called by the main integrator.

3. Configure the code with ``cmake``. The simplest way to trigger *Idefix* configuration is to run ``cmake $IDEFIX_DIR`` in the problem directory, which will genrerate a makefile. More configuration options are described in :ref:`configurationOptions`.

4. Compile the code with ``make``.

5. Create/edit ``idefix.ini``, the input file. This file is read at runtime *only*. It defines most of the properties of the code: resolution, integrator, output, physical modules.

Usually, one starts a new project by copying an example from the test suite of idefix in ``$IDEFIX_DIR/test``. From there, one modifies the setup for the particular problem at hand.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   reference/definitions.hpp
   reference/setup.cpp
   reference/makefile
   reference/idefix.ini
   reference/commandline
   reference/outputs
