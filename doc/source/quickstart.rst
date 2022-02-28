===================
Quickstart tutorial
===================

Download and install *Idefix*
=============================

One first need to download *Idefix* from the public git. Say you want to install *Idefix* in the directory ``~/src/idefix``, you need to run::

    cd ~/src
    git clone https://<yourFavouriteIdefixRepo> idefix
    cd idefix
    git submodule init
    git submodule update

Note that the ``idefix`` directory is automatically created by ``git clone``. Because Kokkos is embedded as a submodule of *Idefix*, the last two lines will download Kokkos for you and put it
at the right place.

The configuration and compilation of *Idefix* relies on the environement variable ``$IDEFIX_DIR``, which should point to the
install directory of *Idefix*. We therefore conclude the installation with::

    export IDEFIX_DIR=~/src/idefix

.. tip::
    The ``IDEFIX_DIR`` environment variable can also be set in the ``.profile`` file of your home directory so that it is automatically set
    when you launch a terminal



Configure and run the SOD tube test problem
===========================================
The test problem are all located in the ``$IDEFIX_DIR/test`` directory of the distribution. To access the Sod shock tube test, one go to::

    cd $IDEFIX_DIR/test/HD/sod

From there, one sees 3 files and a directory:

``definitions.hpp``
    The configuration file. It contains C preprocessor directives for configuration options which require a re-compilation of the code. These are mostly
    the number of vector components, the number of dimensions of the problem, the geometry, etc.

``setup.cpp``
    The physical setup definition. Its role is to define the initial conditions of the problem, but also to provide functions for user-defined
    properties (boundary conditions, diffusivity, potential, etc.)

``idefix.ini``
    The input file. This file is read at runtime *only*. It defines most of the properties of the code: resolution, integrator, output, physical modules.

``python`` directory
    This directory is provided with most of the tests. Its content allows one to check that the code output is consistent with what is expected.

For the time being, the files are already set up for the Sod test problem. The only thing lacking is a ``makefile`` to actually compile the code.
In *Idefix* the makefile is created by `Cmake <https://cmake.org>`_ ,a tool to control code generation on diverse platforms. To configure *Idefix*,
you need Cmake version >=3.16 installed on your machine. For this quickstart, let us configure the code to run on
the cpu in serial (default behaviour). Assuming a ``cmake`` is in the PATH, we simply type in::

    cmake $IDEFIX_DIR

which automatically setups a build tree from the sources in $IDEFIX_DIR using the setup in the current directory.

.. tip::
    If you want to use a specific C++ compiler which is not the default one on your configuration, you can add the ``-D CMAKE_CXX_COMPILER=foo`` option to cmake.

Finally, we compile and run the code::

    make -j 8
    ./idefix

This test being one dimensional, it runs very fast. We can check that the final solution match the prediction of the shock tube problem. To this end, we go to the ``python``
subdirectory and run the test::

    cd python
    python3 ./testidefix.py

If everything goes well, ``testidefix.py`` will load the latest output produced by idefix, display it, compare it with an analytical solution and tell you
whether the error is acceptable or not.


Configure and run the Orszag-Tang test problem
==============================================
The Orszag-Tang problem is a well known test problem for MHD codes. The configuration for this test may be found in::

    cd $IDEFIX_DIR/test/MHD/OrszagTang

As in the Sod test problem, there are 3 files in that directory which completely define the Orszag Tang test problem. We now need to configure the
test with::

    cmake $IDEFIX_DIR

Once the code is configured, it can be ran::

    ./idefix

This test can take much more time that the sod test since it is 2D. In the end, you can visualize the results, which are written as VTK files using
`Paraview <https://www.paraview.org/>`_ of `Visit <https://wci.llnl.gov/simulation/computer-codes/visit>`_. As in the Sod tube test, a ``python/testidefix.py`` script is provided
to check that the last output is consistent with the reference output.

.. tip::
    Given that the Orszag-Tang test can take a long time, you may want to accelerate your computation with a little bit of parallelisation. This can be done with openmp (assuming you have an openmp-compatible compiler)::

        cmake $IDEFIX_DIR -DKokkos_ENABLE_OPENMP=ON
        make -j 8
        export OMP_NUM_THREADS=4
        ./idefix

    or assuming a MPI library is installed on your machine::

        cmake $IDEFIX_DIR -DIdefix_MPI=ON
        make -j 8
        mpirun -np 4 ./idefix

    In both cases, this will run the Orszag-Tang test with 4 threads/processes.
