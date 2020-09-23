===================
Code configuration
===================
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

Problem header file
===================
The problem header file ``definitions.hpp`` contains major parameters which strongly affect the code structure. For this reason, they require re-compilation. A particular effort has been
made to limit their number to a minimum. Each option is a C preprocessor directive. The available options are listed below

``COMPONENTS``
    Sets the number of vector components to the system (for the velocity  and magnetic field). Valid options are 1, 2 or 3.

``DIMENSIONS``
    Sets the number of dimensions of the physical system. Note that ``DIMENSIONS<=COMPONENTS``. It is therefore possible to have more ``COMPONENTS`` than ``DIMENSIONS`` (as is the case
    for axisymmetric systems which have ``COMPONENTS=3`` and ``DIMENSIONS=2``).

``GEOMETRY``
    Sets the global geometry of the coordinate system used by the code. Note that the coordinates used by *Idefix* :math:`(x_1, x_2, x_3)` are fully generic. It is the geometry which
    translates these generic coordinates into physical coordinates

    The available options are:
     + ``CARTESIAN``: a cartesian frame :math:`(x_1,x_2,x_3)=(x,y,z)`
     + ``CYLINDRICAL``: cylindrical frame, which should be used *only* in conjonction with ``DIMENSIONS=2`` (2D axisymmetric setup): :math:`(x_1,x_2,x_3)=(R,z,\phi)`. Note that because the frame is right handed, :math:`\phi` is reversed compared to its usual defitinition. 
     + ``POLAR``: a polar frame: :math:`(x_1,x_2,x_3)=(R,\phi,z)`. This geometry should be used when doing 3D cylindrical problems.
     + ``SPHERICAL``: spherical frame :math:`(x_1,x_2,x_3)=(R,\theta,\phi)` which can be used in 1D, 2D or 3D.


``HAVE_ENERGY``
    Whether or not we solve the total energy equation. If defined, the code uses an ideal equation of state with the adiabatic index :math:`\gamma` set in the input file. Otherwise, the equation of state is isothermal.

``ORDER``
    Spatial order of the reconstruction scheme. As of now, can be 1 (first order, donnor cell reconstruction) or 2 (second order, slope limited Van-leer linear reconstruction)

``EMF_AVERAGE``
    Choose the type of cell corner EMF reconstruction scheme when using MHD. Can be either ``ARITHMETIC``, ``UCT0`` or ``UCT_CONTACT``, following the definition of Gardiner & Stone (2005).

Code configuration
==================

Because the code can be configured for many architectures, it relies on a python configuration script ``$IDEFIX_DIR/configure.py`` to generate the makefile needed. This script accepts
many options to adapt the generated makefile to the architecture on which one wants to run. A complete list of options can be obtained by running ``$IDEFIX_DIR/configure.py -h``. These options are:
``-h, --help``
    Display the help message and exit
``-mhd``
    Enable MHD in the code
``-gpu``
    Enable GPU target architecture (otherwise default target is cpu).
``-arch=xxx``
    Compile for a specific CPU or GPU target. These corresponds to Kokkos target, so user can report to Kokkos documentation to get an up-to-date list of targets. At the time of writing, valid options are 
     + Intel CPUs:    KNC, KNL, SNB, HSW, BDW, SKX
     + NVIDIA GPUs :  Kepler, Kepler30, Kepler32, Kepler35, Kepler37, Maxwell, Maxwell50, Maxwell52, Maxwell53, Pascal60, Pascal61, Volta70, Volta72, Turing75, Ampere80
     + ARM CPUS:      ARMv80, ARMv81, ARMv8-ThunderX, ARMv8-TX2
     + IBM:      BGQ, Power7, Power8, Power9
     + AMD-GPUS: Vega900, Vega906
     + AMD-CPUS: AMDAVX, Zen, Zen2
``-openmp``
    Enable OpenMP parallelisation on supported compilers (not available on GPUs for obvious reasons).

``-mpi``
    Enable MPI (message passing interface) when available. Note that this option can be used in conjonction with -gpu to run idefix simultaneously on several GPUs. This feature requires a CUDA-aware installation of MPI, such as OpenMPI.

Problem Setup
=============


Problem input file
==================


Migrating from PLUTO
====================
