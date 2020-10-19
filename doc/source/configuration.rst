=====================
General configuration
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

Problem header file ``definitions.hpp``
=======================================
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


``ISOTHERMAL``
    If ISOTHERMAL is defined, then the code does not solve the energy equation and instead assume a globally isothermal flow. Otherwise, the flow is assumed to be adiabatic.

``ORDER``
    Spatial order of the reconstruction scheme. As of now, can be 1 (first order, donor cell reconstruction) or 2 (second order, slope limited Van-leer linear reconstruction)

``EMF_AVERAGE``
    Choose the type of cell corner EMF reconstruction scheme when using MHD. Can be either ``ARITHMETIC``, ``UCT0`` or ``UCT_CONTACT``, following the definition of Gardiner & Stone (2005).

Creating the Makefile
=====================

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

Problem Setup ``setup.cpp``
===========================

The setup class
---------------



Function enrollment
-------------------


Problem input file ``idefix.ini``
=================================

The problem input file is by default named ``idefix.ini``. It is possible to start idefix with an other input file using the `-i` command line option.

The problem input file is read when *Idefix* starts. It is splitted into several sections, each section name corresponding to a C++ class in Idefix structure. Inside each section, each line defines an entry, which can have has many parameters as one wish
(note that it requires at least one parameter). The input file
allows for comments, which should start with ``#``.

.. tip::
    Note that you can add arbitray sections and entries in the input file freely. *Idefix* will automatically read and store them on startup. They are then accessible in the code using the
    ``Input::GetReal(..)``, ``Input::GetInt(...)`` and ``Input::GetString(..)`` methods defined in the ``Input`` class.

The ``Grid`` section
--------------------
The grid section defines the grid total dimension. It consists of 3 entries ``X1-grid``, ``X2-grid`` and ``X3-grid``. Each entry defines the repartition of the grid point in the corresponding direction (the grid is always rectilinear).
Each entry defines a series of grid blocks, which can have various spacing. The definition of the grid points is as follows

+-------------+-------------+---------------------+--------------------------+---------------------------------+---------------------------------+----------------------------------------------+-----+---------------------+
|             | Entry name  |   number of blocks  |  begining of first block | number of points in first block | grid spacing in first block     | end of first block/beginning of second block | ... | end of nth block    |
+=============+=============+=====================+==========================+=================================+=================================+==============================================+=====+=====================+
|             | X1/2/3-Grid |  integer number >= 1| floating point           | integer                         | can be u, l, s                  |  floating point                              | ... | floating point      |
+-------------+-------------+---------------------+--------------------------+---------------------------------+---------------------------------+----------------------------------------------+-----+---------------------+
| Example     | X1-Grid     |  1                  |  0.0                     | 64                              |  u                              | 1.0                                          |     |                     |
+-------------+-------------+---------------------+--------------------------+---------------------------------+---------------------------------+----------------------------------------------+-----+---------------------+

In the example above, we define in ``X1`` direction a uniform grid of 64 points starting at ``X1=0.0`` and ending at ``X1=1.0``. The grid spacing can be either uniform (``u``), increasing logarithmically (``l``) or stretched (``s``).

.. warning::
    As of version 0.5, *Idefix* supports only uniform and logarithmically-spaced grids.


The ``TimeIntegrator`` section
------------------------------

This section is used by *Idefix* time integrator class to define the time integrator method and related variables. The entries of this section are as followed


+----------------+--------------------+------------------------------------------------+
|  Entry name    | Parameter type     | Comment                                        |
+================+====================+================================================+
| CFL            | float              | CFL number. Should be < 1 to ensure stability  |
+----------------+--------------------+------------------------------------------------+
| CFL_max_var    | float              | | fraction by which dt is allowed to increase  |
|                |                    | | between two successive timesteps             |
+----------------+--------------------+------------------------------------------------+
| tstop          | float              | time when the code stops                       |
+----------------+--------------------+------------------------------------------------+
| first_dt       | float              | first timestep used by the integrator          |
+----------------+--------------------+------------------------------------------------+
| nstages        | integer            | | number of stages of the integrator. Can be   |
|                |                    | | either 1, 2 or 3. 1=First order Euler method |
|                |                    | | 2, 3= second and third order TVD Runge-Kutta |
+----------------+--------------------+------------------------------------------------+

.. note::
    The ``first_dt`` is necessary since wave speeds are evaluated when Riemann problems are solved, hence the CFL
    condition can only be evaluated after the first timestep.

.. warning::
    As of version 0.4, *Idefix* ignores ``CFL_max_var``, which is by default set to 1.1.


The ``Hydro`` section
---------------------

This section is used by the hydrodynamics class of *Idefix*. It defines the hydrodynamic parameters, and allows one to add some physics. The parameters are as followed:

+----------------+--------------------+----------------------------------------------------------+
|  Entry name    | Parameter type     | Comment                                                  |
+================+====================+==========================================================+
| Solver         | string             | | Type of Riemann Solver. In hydro can be any of         |
|                |                    | | ``tvdlf``, ``hll``, ``hllc`` and ``roe``.              |
|                |                    | | In MHD, can be ``tvdlf``, ``hll``, ``hlld``            |
|                |                    | | and ``roe``                                            |
+----------------+--------------------+----------------------------------------------------------+
| csiso          | float              | | Isothermal sound speed. Only used when                 |
|                |                    | | ISOTHERMAL is defined in ``definition.hpp``            |
+----------------+--------------------+----------------------------------------------------------+
| gamma          | float              | Adiabatic index when ISOTHERMAL is not defined           |
+----------------+--------------------+----------------------------------------------------------+
| Resistivity    | string, float      | | Switches on Ohmic diffusion. String can be             |
|                |                    | | either ``constant`` or ``userdef``.                    |
|                |                    | | When ``constant``, the second parameter is the         |
|                |                    | | Ohmic diffusion coefficient.                           |
|                |                    | | When ``userdef``, the ``Hydro`` class expects a        |
|                |                    | | user-defined diffusivity function to be enrolled with  |
|                |                    | | ``Hydro::EnrollOhmicDiffusivity(DiffusivityFunc)``     |
|                |                    | | In this case, the second parameter is not used.        |
+----------------+--------------------+----------------------------------------------------------+
| Ambipolar      | string, float      | | Switches on ambipolar diffusion. String can be         |
|                |                    | | either ``constant`` or ``userdef``.                    |
|                |                    | | When ``constant``, the second parameter is the         |
|                |                    | | ambipolar diffusion coefficient.                       |
|                |                    | | When ``userdef``, the ``Hydro`` class expects a        |
|                |                    | | user-defined diffusivity function to be enrolled with  |
|                |                    | | ``Hydro::EnrollAmbipolarDiffusivity(DiffusivityFunc)`` |
|                |                    | | In this case, the second parameter is not used.        |
+----------------+--------------------+----------------------------------------------------------+
| Hall           | string, float      | | Switches on Hall effect. String can be                 |
|                |                    | | either ``constant`` or ``userdef``.                    |
|                |                    | | When ``constant``, the second parameter is the         |
|                |                    | | Hall diffusion coefficient.                            |
|                |                    | | When ``userdef``, the ``Hydro`` class expects a        |
|                |                    | | user-defined diffusivity function to be enrolled with  |
|                |                    | | ``Hydro::EnrollHallDiffusivity(DiffusivityFunc)``      |
|                |                    | | In this case, the second parameter is not used.        |
+----------------+--------------------+----------------------------------------------------------+
| GravPotential  | string             | | Switches on an external gravitational potential.       |
|                |                    | | Only ``userdef`` is allowed.                           |
|                |                    | | When ``userdef is set, the ``Hydro`` class expects     |
|                |                    | | a user-defined potential function to be enrolled with  |
|                |                    | | ``Hydro::EnrollGravPotential(GravPotentialFunc)``      |   
+----------------+--------------------+----------------------------------------------------------+
| Rotation       | float,float,float  | | Add rotation with rhe rotation vector components given |
|                |                    | | as parameters. This entry only adds Coriolis force.    |
+----------------+--------------------+----------------------------------------------------------+
| ShearingBox    | float              | | Enable shearing box source terms.  The entry parameter |
|                |                    | | corresponds to the shear rate                          |
|                |                    |  :math:`dv_{x2}/d x_1`.                                  |
|                |                    | | Note that this is not sufficient to fully define a     |
|                |                    | | shearing box: boundary conditions are also required.   | 
+----------------+--------------------+----------------------------------------------------------+



.. note::
    The Hall effect is implemented directly in the HLL Riemann solver following Lesur, Kunz & Fromang (2014)
    and adding the whistler speed only to the magnetic flux function, following Marchand et al. (2019).
    For these reasons, Hall can only be used in conjonction with the HLL Riemann solver. In addition, only
    the arithmetic Emf reconstruction scheme has been shown to work systematically with Hall, and is therefore
    strongly recommended for production runs.   

The ``Boundary`` section
------------------------

This section describes the boundary conditions of the code. 


Command line options
====================


Migrating from PLUTO
====================
