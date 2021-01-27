Problem input file ``idefix.ini``
=================================

The problem input file is by default named ``idefix.ini``. It is possible to start idefix with an other input file using the `-i` command line option.

The problem input file is read when *Idefix* starts. It is splitted into several sections, each section name corresponding to a C++ class in Idefix structure. Inside each section, each line defines an entry, which can have as many parameters as one wishes
(note that it requires at least one parameter). The input file
allows for comments, which should start with ``#``.

.. tip::
    Note that you can add arbitray sections and entries in the input file freely. *Idefix* will automatically read and store them on startup. They are then accessible in the code using the
    ``Input::GetReal(..)``, ``Input::GetInt(...)`` and ``Input::GetString(..)`` methods defined in the ``Input`` class (see :ref:`inputClass`)

``Grid`` section
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

In the example above, we define in ``X1`` direction a uniform grid (``u``) of 64 points starting at ``X1=0.0`` and ending at ``X1=1.0``.

The grid spacing can be one of the following:

* Uniform spacing (``u``): a block with constant spacing is constructed. The grid spacing :math:`\Delta x` is defined as :math:`\Delta x=\frac{x_\mathrm{end}-x_\mathrm{start}}{N}`

* Logarthmic spacing  (``l``): the block spacing is proportional to the coordinate :math:`\Delta x\propto x`. More formally, the cell boundaries are defined as  :math:`x_{i-1/2}=x_\mathrm{start}\alpha^{i/N}` where  :math:`\alpha=\frac{x_\mathrm{end}+|x_\mathrm{start}|-x_\mathrm{start}}{|x_\mathrm{start}|}`. The grid spacing is then defined through :math:`\Delta x_i=x_{i+1/2}-x_{i-1/2}`.

* Stretched spacing (``s+`` or ``s-``): the block spacing is defined to follow a geometrical series, starting from the reference spacing :math:`\Delta x_0` taken from the previous (``-``) or next (``+``) uniform block. Mathematically, the grid spacing is defined as :math:`\Delta x_i=\Delta x_0 r^{i-1}` for ``s+`` and  :math:`\Delta x_i=\Delta x_0 r^{N-i}` for ``s-``. The streching ratio :math:`r` is itself defined implicitly as :math:`\frac{r-r^{N+1}}{1-r}=\frac{x_\mathrm{end}-x_\mathrm{start}}{\Delta x_0}`


.. tip::
  Note that the stretched block requires at least one uniform block on one of it side to define the reference spacing :math:`\Delta x_0`

``TimeIntegrator`` section
------------------------------

This section is used by *Idefix* time integrator class to define the time integrator method and related variables. The entries of this section are as followed


+----------------+--------------------+-----------------------------------------------------------------------------------------------------------+
|  Entry name    | Parameter type     | Comment                                                                                                   |
+================+====================+===========================================================================================================+
| CFL            | float              | CFL number. Should be < 1 to ensure stability                                                             |
+----------------+--------------------+-----------------------------------------------------------------------------------------------------------+
| CFL_max_var    | float              | fraction by which :math:`dt` is allowed to increase between two  successive timesteps                     |
+----------------+--------------------+-----------------------------------------------------------------------------------------------------------+
| tstop          | float              | time when the code stops                                                                                  |
+----------------+--------------------+-----------------------------------------------------------------------------------------------------------+
| first_dt       | float              | first timestep used by the integrator                                                                     |
+----------------+--------------------+-----------------------------------------------------------------------------------------------------------+
| nstages        | integer            | | number of stages of the integrator. Can be  either 1, 2 or 3. 1=First order Euler method,               |
|                |                    | | 2, 3 = second and third order  TVD Runge-Kutta                                                          |
+----------------+--------------------+-----------------------------------------------------------------------------------------------------------+

.. note::
    The ``first_dt`` is necessary since wave speeds are evaluated when Riemann problems are solved, hence the CFL
    condition can only be evaluated after the first timestep.


``Hydro`` section
---------------------

This section is used by the hydrodynamics class of *Idefix*. It defines the hydrodynamic parameters, and allows one to add some physics. The parameters are as followed:

+----------------+-------------------------+---------------------------------------------------------------------------------------------+
|  Entry name    | Parameter type          | Comment                                                                                     |
+================+=========================+=============================================================================================+
| Solver         | string                  | | Type of Riemann Solver. In hydro can be any of ``tvdlf``, ``hll``, ``hllc`` and ``roe``.  |
|                |                         | | In MHD, can be ``tvdlf``, ``hll``, ``hlld`` and ``roe``                                   |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| csiso          | string, float           | | Isothermal sound speed. Only used when ISOTHERMAL is defined in ``definitions.hpp``.      |
|                |                         | | When ``constant``, the second parameter is the spatially constant sound speed.            |
|                |                         | | When ``userdef``, the ``Hydro`` class expects a user-defined sound speed function         |
|                |                         | | to be enrolled with   ``EnrollIsoSoundSpeed(IsoSoundSpeedFunc)``                          |
|                |                         | | (see :ref:`functionEnrollment`). In this case, the second parameter is not used.          |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| gamma          | float                   | Adiabatic index when ISOTHERMAL is not defined                                              |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| Resistivity    | string, float           | | Switches on Ohmic diffusion. String can be  either ``constant`` or ``userdef``.           |
|                |                         | | When ``constant``, the second parameter is the  Ohmic diffusion coefficient.              |
|                |                         | | When ``userdef``, the ``Hydro`` class expects a user-defined diffusivity function         |
|                |                         | | to be enrolled with   ``Hydro::EnrollOhmicDiffusivity(DiffusivityFunc)``                  |
|                |                         | | (see :ref:`functionEnrollment`). In this case, the second parameter is not used.          |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| Ambipolar      | string, float           | | Switches on ambipolar diffusion. String can be  either ``constant`` or ``userdef``.       |
|                |                         | | When ``constant``, the second parameter is the ambipolar diffusion coefficient.           |
|                |                         | | When ``userdef``, the ``Hydro`` class expects a user-defined diffusivity function         |
|                |                         | | to be enrolled with   ``Hydro::EnrollAmbipolarDiffusivity(DiffusivityFunc)``              |
|                |                         | | (see :ref:`functionEnrollment`). In this case, the second parameter is not used.          |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| Hall           | string, float           | | Switches on Hall effect. String can be  either ``constant`` or ``userdef``.               |
|                |                         | | When ``constant``, the second parameter is the  Hall diffusion coefficient.               |
|                |                         | | When ``userdef``, the ``Hydro`` class expects a user-defined diffusivity function         |
|                |                         | | to be enrolled with   ``Hydro::EnrollHallDiffusivity(DiffusivityFunc)``                   |
|                |                         | | (see :ref:`functionEnrollment`). In this case, the second parameter is not used.          |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| Viscosity      | string, float, float    | | Switches on viscous diffusion. String can be either ``constant`` or ``userdef``           |
|                |                         | | When ``constant``, the second parameter is the flow viscosity and the third               |
|                |                         | | parameter is the second (or compressive) viscosity (which is optionnal).                  |
|                |                         | | When ``userdef``, the ``Hydro.Viscosity`` class expects a user-defined viscosity function |
|                |                         | | to be enrolled with   ``Hydro.Viscosity::EnrollViscousDiffusivity(DiffusivityFunc)``      |
|                |                         | | (see :ref:`functionEnrollment`). In this case, the second and third parameters            |
|                |                         | | are not used.                                                                             |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| GravPotential  | string                  | | Switches on an external gravitational potential. Only ``userdef`` is allowed.             |
|                |                         | | When ``userdef is set, the ``Hydro`` class expects  a user-defined potential function     |
|                |                         | | to be enrolled with  ``Hydro::EnrollGravPotential(GravPotentialFunc)``                    |
|                |                         | | (see :ref:`functionEnrollment`)                                                           |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| Rotation       | float,float,float       | | Add rotation with rhe rotation vector components given as parameters.                     |
|                |                         | | Note that this entry only adds Coriolis force.                                            |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| ShearingBox    | float                   | | Enable shearing box source terms.  The entry parameter corresponds to the shear rate      |
|                |                         | | :math:`dv_{x2}/d x_1`.                                                                    |
|                |                         | | Note that this is not sufficient to fully define a shearing box: boundary conditions      |
|                |                         | | are also required.                                                                        |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+



.. note::
    The Hall effect is implemented directly in the HLL Riemann solver following Lesur, Kunz & Fromang (2014)
    and adding the whistler speed only to the magnetic flux function, following Marchand et al. (2019).
    For these reasons, Hall can only be used in conjonction with the HLL Riemann solver. In addition, only
    the arithmetic Emf reconstruction scheme has been shown to work systematically with Hall, and is therefore
    strongly recommended for production runs.

``Boundary`` section
------------------------

This section describes the boundary conditions used by the code. There are 6 entries
which need to be defined: ``X1-beg``, ``X2-beg``, ``X3-beg`` for the left boundaries in the direction X1, X2, X3,
and ``X1-end``, ``X2-end``, ``X3-end`` for the right boundaries. Each boundary can be assigned the following types of conditions

+----------------+------------------------------------------------------------------------------------------------------------------+
| Boundary type  | Comment                                                                                                          |
+================+==================================================================================================================+
| outflow        | | zero gradient on the density, pressure, tangential velocity and magnetic field. The normal velocity is set to  |
|                | | zero gradient when it is flowing outwards otherwise it is set to 0.                                            |
+----------------+------------------------------------------------------------------------------------------------------------------+
| periodic       |  Periodic boundary conditions. Each field is copied between beg and end sides of the boundary.                   |
+----------------+------------------------------------------------------------------------------------------------------------------+
| reflective     | The normal component of the velocity is systematically reversed. Otherwise identical to ``outflow``.             |
+----------------+------------------------------------------------------------------------------------------------------------------+
| shearingbox    | Shearing-box boudary conditions.                                                                                 |
+----------------+------------------------------------------------------------------------------------------------------------------+
| userdef        | | User-defined boundary conditions. The boundary condition function should be enrolled in the setup constructor  |
|                | | (see :ref:`userdefBoundaries`)                                                                                 |
+----------------+------------------------------------------------------------------------------------------------------------------+

.. warning::
    As of version 0.5, *Idefix* shearing-box assumes axisymmetry, i.e. NX2=1.


.. _outputSection:

``Output`` section
----------------------

This section describes the outputs *Idefix* produces. For more details about each output type, have a look at :ref:`output`.

+----------------+-------------------------+--------------------------------------------------------------------------------------------------+
|  Entry name    | Parameter type          | Comment                                                                                          |
+================+=========================+==================================================================================================+
| log            | integer                 | | Time interval between log outputs, in code steps (default 100).                                |
+----------------+-------------------------+--------------------------------------------------------------------------------------------------+
| dmp            | float                   | | Time interval between dump outputs, in code units.                                             |
|                |                         | | If negative, periodic dump outputs are disabled.                                               |
+----------------+-------------------------+--------------------------------------------------------------------------------------------------+
| vtk            | float                   | | Time interval between vtk outputs, in code units.                                              |
|                |                         | | If negative, periodic vtk outputs are disabled.                                                |
+----------------+-------------------------+--------------------------------------------------------------------------------------------------+
| analysis       | float                   | | Time interval between analysis outputs, in code units.                                         |
|                |                         | | If negative, periodic analysis outputs are disabled.                                           |
|                |                         | | When this entry is set, *Idefix* expects a user-defined analysis function to be                |
|                |                         | | enrolled with  ``Output::EnrollAnalysis(AnalysisFunc)`` (see :ref:`functionEnrollment`).       |
+----------------+-------------------------+--------------------------------------------------------------------------------------------------+
| uservar        | string series           | | List the name of the user-defined variables the user wants to define.                          |
|                |                         | | When this list is present in the input file, *Idefix* expects a user-defined                   |
|                |                         | | function to be enrolled with ``Output::EnrollUserDefVariables(UserDefVariablesFunc)``          |
|                |                         | | (see :ref:`functionEnrollment`). The user-defined variables defined by this function           |
|                |                         | | are then written as new variables in vtk  outputs.                                             |
+----------------+-------------------------+--------------------------------------------------------------------------------------------------+

.. note::
    Even if dumps are not mentionned in your input file (and are therefore disabled), dump files are still produced when *Idefix* captures a signal
    (see :ref:`signalHandling`).