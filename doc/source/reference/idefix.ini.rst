Problem input file ``idefix.ini``
=================================

The problem input file is by default named ``idefix.ini``. It is possible to start idefix with an other input file using the `-i` command line option.

The problem input file is read when *Idefix* starts. It is split into several sections, each section name corresponding to a C++ class in Idefix structure. Inside each section, each line defines an entry, which can have as many parameters as one wishes
(note that it requires at least one parameter). The input file
allows for comments, which should start with ``#``.

.. tip::
    You can add arbitray sections and entries in the input file freely. *Idefix* will automatically read and store them on startup. They are then accessible in the code using the
    ``Input::Get<T>(..)`` and ``Input::GetOrSet<T>`` template methods defined in the ``Input`` class (see :ref:`inputClass`).
    To avoid any name collisions with future versions of Idefix, we recommend setting setup-specific parameters in a ``[Setup]`` section.

.. _gridSection:

``Grid`` section
--------------------
The grid section defines the grid total dimension. It consists of 3 entries ``X1-grid``, ``X2-grid`` and ``X3-grid``. Each entry defines the repartition of the grid points in the corresponding direction (the grid is always rectilinear).
Each entry defines a series of grid blocks which are concatenated along the direction. Each block in a direction can have a different spacing rule (uniform, log or stretched). The definition of the Grid entries is as follows

+----------------------------+-------------------------+------------------------------+
|                            |  Allowed value          |    Example                   |
+============================+=========================+==============================+
| Entry name                 | X1/2/3-Grid             | X1-Grid                      |
+----------------------------+-------------------------+------------------------------+
| # of grid blocks           | integer >= 1            | 1                            |
+----------------------------+-------------------------+------------------------------+
| start of 1st block         | floating point          | 0.0                          |
+----------------------------+-------------------------+------------------------------+
| # of points in 1st block   | integet >= 1            | 64                           |
+----------------------------+-------------------------+------------------------------+
| spacing in first block     | "u", "l", "s+", "s-"    | u                            |
+----------------------------+-------------------------+------------------------------+
| | end of 1st block         | floating point          | 1.0                          |
| | /start of 2nd block      |                         |                              |
+----------------------------+-------------------------+------------------------------+
| ...repeat for each block...| ...                     | ...                          |
+----------------------------+-------------------------+------------------------------+

In the example above, we define in ``X1`` direction a uniform grid (``u``) of 64 points starting at ``X1=0.0`` and ending at ``X1=1.0``.
This would be written in the input file as:

.. code-block::

  [Grid]
  X1-Grid        1     0.0   64     u     1.0


The grid spacing can be one of the following:

* Uniform spacing (``u``): a block with constant spacing is constructed. The grid spacing :math:`\Delta x` is defined as :math:`\Delta x=\frac{x_\mathrm{end}-x_\mathrm{start}}{N}`

* Logarthmic spacing  (``l``): the block spacing is proportional to the coordinate :math:`\Delta x\propto x`. More formally, the cell boundaries are defined as  :math:`x_{i-1/2}=x_\mathrm{start}\alpha^{i/N}` where  :math:`\alpha=\frac{x_\mathrm{end}+|x_\mathrm{start}|-x_\mathrm{start}}{|x_\mathrm{start}|}`. The grid spacing is then defined through :math:`\Delta x_i=x_{i+1/2}-x_{i-1/2}`.

* Stretched spacing (``s+`` or ``s-``): the block spacing is defined to follow a geometrical series, starting from the reference spacing :math:`\Delta x_0` taken from the previous (``-``) or next (``+``) uniform block. Mathematically, the grid spacing is defined as :math:`\Delta x_i=\Delta x_0 r^{i-1}` for ``s+`` and  :math:`\Delta x_i=\Delta x_0 r^{N-i}` for ``s-``. The streching ratio :math:`r` is itself defined implicitly as :math:`\frac{r-r^{N+1}}{1-r}=\frac{x_\mathrm{end}-x_\mathrm{start}}{\Delta x_0}`


.. note::
  Note that the stretched block requires at least one uniform block on one of it side to define the reference spacing :math:`\Delta x_0`

.. tip::
  It is also possible to change the grid spacing to increase the integration timestep with the ``coarsening`` entry, which enables grid coarsening
  (see :ref:`gridCoarseningModule`)

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
| first_dt       | float              | first timestep used by the integrator. If not set, Idefix use by default 1e-10 (very conservative)        |
+----------------+--------------------+-----------------------------------------------------------------------------------------------------------+
| fixed_dt       | float              | | when set, *Idefix* uses a fixed time step instead of the dt computed from the CFL condition.            |
|                |                    | | In this case, the CFL parameters and ``first_dt`` are ignored.                                          |
+----------------+--------------------+-----------------------------------------------------------------------------------------------------------+
| max_runtime    | float              | | when set, *Idefix* aborts the calculation when it has run for `max_runtime` hours (wall clock time).    |
|                |                    | | In this case, a restart dump is automatically written when the code stops.                              |
+----------------+--------------------+-----------------------------------------------------------------------------------------------------------+
| nstages        | integer            | | number of stages of the integrator. Can be  either 1, 2 or 3. 1=First order Euler method,               |
|                |                    | | 2, 3 = second and third order  TVD Runge-Kutta                                                          |
+----------------+--------------------+-----------------------------------------------------------------------------------------------------------+
| check_nan      | integer            | | number of time integration cycles between each Nan verification. Default is 100.                        |
|                |                    | | Note that Nan checks are slow on GPUs, and low values of ``check_nan`` are not recommended.             |
+----------------+--------------------+-----------------------------------------------------------------------------------------------------------+
| maxdivB        | float              |  Maximum divB tolerated. Default is 1e-6 in double precision and 1e-2 in single precision.                |
+----------------+--------------------+-----------------------------------------------------------------------------------------------------------+

.. note::
    The ``first_dt`` is recommended since wave speeds are evaluated when Riemann problems are solved, hence the CFL
    condition can only be evaluated after the first timestep.


``Hydro`` section
---------------------

This section is used by the hydrodynamics class of *Idefix*. It defines the hydrodynamic parameters, and allows one to add some physics. The parameters are as followed:

+----------------+-------------------------+---------------------------------------------------------------------------------------------+
|  Entry name    | Parameter type          | Comment                                                                                     |
+================+=========================+=============================================================================================+
| solver         | string                  | | Type of Riemann Solver. In hydro can be any of ``tvdlf``, ``hll``, ``hllc`` and ``roe``.  |
|                |                         | | In MHD, can be ``tvdlf``, ``hll``, ``hlld`` and ``roe``                                   |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| emf            | string                  | | Averaging scheme for the electromotive force (only used with MHD). The options            |
|                |                         | | follows Gardiner & Stone JCP, 2005 (GS05).                                                |
|                |                         | | ``arithmetic``: simple arithmetic average of the face-centered emfs (eq. 33 in GS05)      |
|                |                         | | ``uct0``: Upwind constraint transport (UCT) with 0 wave speed (eq. 39 in GS05)            |
|                |                         | | ``uct_contact``: UCT with contact wave upwinding (eq. 50 in GS05)                         |
|                |                         | | ``uct_hll``: UCT with 2D Riemann solver using the HLL approximation. Follows Londrillo    |
|                |                         | |  & del Zanna JCP (2004).                                                                  |
|                |                         | | ``uct_hlld``: UCT with 2D Riemann solver using the HLLD approximation. Follows Londrillo  |
|                |                         | |  & del Zanna JCP (2004).                                                                  |
|                |                         | |  If no averaging scheme is selected in the input file, *Idefix* uses ``uct_contact``.     |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| csiso          | string, (float)         | | Isothermal sound speed. Only used when ISOTHERMAL is defined in ``definitions.hpp``.      |
|                |                         | | When ``constant``, the second parameter is the spatially constant sound speed.            |
|                |                         | | When ``userdef``, the ``Hydro`` class expects a user-defined sound speed function         |
|                |                         | | to be enrolled with   ``EnrollIsoSoundSpeed(IsoSoundSpeedFunc)``                          |
|                |                         | | (see :ref:`functionEnrollment`). In this case, the second parameter is not used.          |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| gamma          | float                   | | Adiabatic index when ISOTHERMAL is not defined. Default to 5/3 if not set.                |
|                |                         | | NB: this parameter is used only by the default equation of state implemented in *Idefix*  |
|                |                         | | Custom equation of states (:ref:`eosModule`) ignore this parameter                        |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| tracer         | integer                 | Number of passive tracers associated to the fluid. Default to 0 if not set.                 |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| resistivity    | string, string, (float) | | Switches on Ohmic diffusion.                                                              |
|                |                         | | The first parameter can be ``explicit`` or ``rkl``. When ``explicit``, diffusion is       |
|                |                         | | integrated in the main integration loop with the usual cfl restriction.  If ``rkl``,      |
|                |                         | | diffusion  is integrated using the Runge-Kutta Legendre scheme.                           |
|                |                         | | The second String can be  either ``constant`` or ``userdef``.                             |
|                |                         | | When ``constant``, the second parameter is the  Ohmic diffusion coefficient.              |
|                |                         | | When ``userdef``, the ``Hydro`` class expects a user-defined diffusivity function         |
|                |                         | | to be enrolled with   ``Hydro::EnrollOhmicDiffusivity(DiffusivityFunc)``                  |
|                |                         | | (see :ref:`functionEnrollment`). In this case, the third  parameter is not used.          |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| ambipolar      | string, string, (float) | | Switches on ambipolar diffusion.                                                          |
|                |                         | | The first parameter can be ``explicit`` or ``rkl``. When ``explicit``, diffusion is       |
|                |                         | | integrated in the main integration loop with the usual cfl restriction.  If ``rkl``,      |
|                |                         | | diffusion  is integrated using the Runge-Kutta Legendre scheme.                           |
|                |                         | | The second String can be  either ``constant`` or ``userdef``.                             |
|                |                         | | When ``constant``, the second parameter is the ambipolar diffusion coefficient.           |
|                |                         | | When ``userdef``, the ``Hydro`` class expects a user-defined diffusivity function         |
|                |                         | | to be enrolled with   ``Hydro::EnrollAmbipolarDiffusivity(DiffusivityFunc)``              |
|                |                         | | (see :ref:`functionEnrollment`). In this case, the third parameter is not used.           |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| hall           | string, string, (float) | | Switches on Hall effect.                                                                  |
|                |                         | | The first parameter can only be ``explicit``.                                             |
|                |                         | | The second String can be  either ``constant`` or ``userdef``.                             |
|                |                         | | When ``constant``, the third parameter is the  Hall diffusion coefficient.                |
|                |                         | | When ``userdef``, the ``Hydro`` class expects a user-defined diffusivity function         |
|                |                         | | to be enrolled with   ``Hydro::EnrollHallDiffusivity(DiffusivityFunc)``                   |
|                |                         | | (see :ref:`functionEnrollment`). In this case, the third parameter is not used.           |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| viscosity      | string, string,         | | Switches on viscous diffusion.                                                            |
|                | float, (float)          | | The first parameter can be ``explicit`` or ``rkl``. When ``explicit``, diffusion is       |
|                |                         | | integrated in the main integration loop with the usual cfl restriction.  If ``rkl``,      |
|                |                         | | diffusion  is integrated using the Runge-Kutta Legendre scheme.                           |
|                |                         | | The second parameter can be  either ``constant`` or ``userdef``.                          |
|                |                         | | When ``constant``, the third parameter is the flow viscosity and the fourth               |
|                |                         | | parameter is the second (or compressive) viscosity (which is optionnal).                  |
|                |                         | | When ``userdef``, the ``Hydro.Viscosity`` class expects a user-defined viscosity function |
|                |                         | | to be enrolled with   ``Hydro.Viscosity::EnrollViscousDiffusivity(DiffusivityFunc)``      |
|                |                         | | (see :ref:`functionEnrollment`). In this case, the third and fourth parameters            |
|                |                         | | are not used.                                                                             |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| TDiffusion     | string, string,         | | Switches on isotropic thermal diffusion.                                                  |
|                | float                   | | The first parameter can be ``explicit`` or ``rkl``. When ``explicit``, diffusion is       |
|                |                         | | integrated in the main integration loop with the usual cfl restriction.  If ``rkl``,      |
|                |                         | | diffusion  is integrated using the Runge-Kutta Legendre scheme.                           |
|                |                         | | The second parameter can be  either ``constant`` or ``userdef``.                          |
|                |                         | | When ``constant``, the third parameter is the (constant) thermal diffusivity.             |
|                |                         | | When ``userdef``, the ``Hydro.ThermalDiffusivity`` class expects a user-defined thermal   |
|                |                         | | diffusivity function to be enrolled with                                                  |
|                |                         | | ``Hydro.thermalDiffusion::EnrollThermalDiffusivity(DiffusivityFunc)`` .                   |
|                |                         | | (see :ref:`functionEnrollment`) In this case, the third parameter is not used.            |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| rotation       | float                   | | Add rotation with the z rotation speed given as parameter.                                |
|                |                         | | Note that this entry only adds Coriolis force in Cartesian geometry.                      |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| shearingBox    | float                   | | Enable shearing box source terms.  The entry parameter corresponds to the shear rate      |
|                |                         | | :math:`dv_{x2}/d x_1`.                                                                    |
|                |                         | | Note that this is not sufficient to fully define a shearing box: boundary conditions      |
|                |                         | | are also required.                                                                        |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| shockFlattening| float                   | | Enable shock flattening.  When enabled, the reconstruction scheme reverts to minmod       |
|                |                         | | limiter when strong shocks are detected. The entry parameter is the threshold above which |
|                |                         | | a shock is considered "strong", in units of :math:`|\nabla P /P|`. A low value hence tends|
|                |                         | | to increase the code numerical diffusivity. Typical values are 1 to 10.                   |
|                |                         | | Note that it is possible to enroll a user-defined function to flag specific cells for     |
|                |                         | | shock flattening, in addition to the default flag. This user function can be enrolled     |
|                |                         | | with ``Hydro.shockFlattening.EnrollUserShockFlag(UserShockFunc)`` .                       |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+


.. note::
    The Hall effect is implemented directly in the HLL Riemann solver following Lesur, Kunz & Fromang (2014)
    and adding the whistler speed only to the magnetic flux function, following Marchand et al. (2019).
    For these reasons, Hall can only be used in conjonction with the HLL Riemann solver. In addition, only
    the arithmetic Emf reconstruction scheme has been shown to work systematically with Hall, and is therefore
    strongly recommended for production runs.

.. _fargoSection:

``Fargo`` section
------------------

This section enables the orbital advection algorithm provided in *Idefix*. More information may be found in :ref:`fargoModule`

+----------------+-------------------------+---------------------------------------------------------------------------------------------+
|  Entry name    | Parameter type          | Comment                                                                                     |
+================+=========================+=============================================================================================+
| velocity       | string                  | | Defines orbital advection (Fargo-like) velocity to speed up integration when a strong     |
|                |                         | | azimuthal motion is present (as in a thin disk).  The ``velocity`` can be either          |
|                |                         | | `shearingbox` or `userdef`.                                                               |
|                |                         | | When `shearingbox`, the fargo module uses the linear shear computed by the shearing box   |
|                |                         | | module as the input velocity function.                                                    |
|                |                         | | When `userdef` is set, the fargo module expects a user-defined  velocity function to      |
|                |                         | | be enrolled via Fargo::EnrollVelocity(FargoVelocityFunc)                                  |
|                |                         | |                                                                                           |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| maxShift       | integer                 | | optional: when using MPI with a domain decomposition in the azimuthal direction, this sets|
|                |                         | | the maximum number of cells Fargo is allowed to shift the domain at each time step.       |
|                |                         | | Default: 10                                                                               |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+

.. _gravitySection:

``Gravity`` section
--------------------

This section enables gravity in the form of a gravitational potential and/or an acceleration vector. The gravitational potential used by the code
reads

:math:`\psi=-G_c M_{\rm central}/R+\psi_{SG}+\psi_{\rm userdef}`

where :math:`G_c` is the gravitational constant, :math:`M_{\rm central}` is the mass of a central object, :math:`\psi_{SG}` is the
self-gravitational potential and :math:`\psi_{\rm userdef}` is a user-defined potential. Each term can be enabled individually in the gravity
section as followed:

+----------------+-------------------------+---------------------------------------------------------------------------------------------+
|  Entry name    | Parameter type          | Comment                                                                                     |
+================+=========================+=============================================================================================+
| potential      | string, [string...]     | | Switches on an external gravitational potential. Each parameter adds a potential to the   |
|                |                         | | total potential used by *Idefix*.                                                         |
|                |                         | |                                                                                           |
|                |                         | | * ``userdef`` allows the user to give *Idefix* a user-defined potential function. In this |
|                |                         | | case, ``Gravity`` class expects a user-defined potential function to be enrolled with     |
|                |                         | | ``Gavity::EnrollPotential(GravPotentialFunc)``  (see :ref:`functionEnrollment`)           |
|                |                         | | * ``central`` allows the user to automatically add the potential of a central point mass. |
|                |                         | | In this case, the central mass is assumed to be 1 in code units. This can be modified     |
|                |                         | | using the Mcentral parameter, or using the ``Gravity::SetCentralMass(real)`` method.      |
|                |                         | | * ``selfgravity`` enables the potential computed from solving Poisson equation with the   |
|                |                         | | density distribution (see :ref:`selfGravitySection` and :ref:`selfGravityModule`).        |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| Mcentral       | real                    | | Mass of the central object when a central potential is enabled (see above). Default is 1. |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| gravCst        | real                    | | Set the value of the gravitational constant :math:`G_c` used by the central               |
|                |                         | | mass potential and self-gravitational potential (when enabled) ). Default is 1.           |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| bodyForce      | string                  | | Adds an acceleration vector to each cell of the domain. The only value allowed            |
|                |                         | | is ``userdef``. The ``Gravity`` class then expects a user-defined bodyforce function to   |
|                |                         | | be enrolled via ``Gavity::EnrollBodyForce(BodyForceFunc)`` (see :ref:`functionEnrollment`)|
|                |                         | | See the shearing box tests for examples of using bodyForce.                               |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| skip           | int                     | | Set the number of integration cycles between each computation of the gravity potential.   |
|                |                         | | Default is 1 (i.e. gravity is computed at every cycle).                                   |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+



.. _selfGravitySection:

``SelfGravity`` section
-----------------------

This section describes the method used to compute the self-gravitating potential :math:`\psi_{SG}`. More details on the algorithm may be found in the dedicated
:ref:`selfGravityModule` documentation. For this module to be used, self-gravity must be enabled as a source
of gravitational potential in the ``Gravity`` section (see :ref:`gravitySection` above).

+----------------+-------------------------+---------------------------------------------------------------------------------------------+
|  Entry name    | Parameter type          | Comment                                                                                     |
+================+=========================+=============================================================================================+
| solver         | string                  | | Specifies which solver should be used. Can be ``Jacobi``, ``BICGSTAB`` or ``PBICGSTAB``   |
|                |                         | | for the left preconditionned BICGSTAB solve.                                              |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| targetError    | real                    | | Set the error allowed in the residual :math:`r=\Delta\psi_G/(4\pi G_c)-\rho`. The error   |
|                |                         | | computation is based on a L2 norm. Default is 1e-2.                                       |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| maxIter        | int                     | | Set the maximum number of iterations allowed to the solver to reach convergence. Default  |
|                |                         | | is 1000.                                                                                  |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| boundary-Xn-dir| string                  | | Boundary condition applied to the potential field computed by self-gravity                |
|                |                         | | ``n`` can be 1, 2 or 3 and is the direction for the boundary condition. ``dir`` can be    |
|                |                         | | ``beg`` or ``end`` and indicates the side of the boundary.                                |
|                |                         | | The boundary conditions allowed by the self-gravity solver are described in               |
|                |                         | | :ref:`selfGravityModule`                                                                  |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| skip           | int                     | | Set the number of integration cycles between each computation of self-gravity potential.  |
|                |                         | | Default is 1 (i.e. self-gravity is computed at every cycle).                              |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+



``RKL`` section
------------------

This section controls the Runge-Kutta-Legendre integration module. RKL is automatically enabled when parabolic terms use the `rkl` option. Otherwise,
this block is simply ignored.

+----------------+--------------------+-----------------------------------------------------------------------------------------------------------+
|  Entry name    | Parameter type     | Comment                                                                                                   |
+================+====================+===========================================================================================================+
| cfl            | float              | CFL number for the RKL sub-step. Should be <0.5 for stability. Set by default to 0.5 if not provided      |
+----------------+--------------------+-----------------------------------------------------------------------------------------------------------+
| rmax_par       | float              | Maximum ratio between the hyperbolic timestep and the parabolic (RKL) timestep. Set to 100.0 by default.  |
+----------------+--------------------+-----------------------------------------------------------------------------------------------------------+
| check_nan      | bool               | Whether RKL should check the solution when running. This option affects performances. Default false.      |
+----------------+--------------------+-----------------------------------------------------------------------------------------------------------+

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
| axis           | | Axis Boundary conditions. Useful if one wants to include the axis in spherical geometry in the computational   |
|                | | domain. This condition explicitely requires X2 to go from 0 to :math:`\pi` but can be used for domains         |
|                | | extending over a fraction of a full circle in X3 (i.e :math:`2\pi/n` where :math:`n` is an integer). When the  |
|                | | X3 domain spans :math:`2\pi` and MPI is used, the number of processes along the X3 direction should be one or  |
|                | | even (in this last case, additional communications are required which may impact performances).                |
+----------------+------------------------------------------------------------------------------------------------------------------+
| userdef        | | User-defined boundary conditions. The boundary condition function should be enrolled in the setup constructor  |
|                | | (see :ref:`userdefBoundaries`)                                                                                 |
+----------------+------------------------------------------------------------------------------------------------------------------+


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
| dmp_dir        | string                  | | directory for dump file outputs. Default to "./"                                               |
|                |                         | | The directory is automatically created if it does not exist.                                   |
+----------------+-------------------------+--------------------------------------------------------------------------------------------------+
| vtk            | float                   | | Time interval between vtk outputs, in code units.                                              |
|                |                         | | If negative, periodic vtk outputs are disabled.                                                |
+----------------+-------------------------+--------------------------------------------------------------------------------------------------+
| vtk_dir        | string                  | | directory for vtk file outputs. Default to "./"                                                |
|                |                         | | The directory is automatically created if it does not exist.                                   |
+----------------+-------------------------+--------------------------------------------------------------------------------------------------+
| vtk_sliceN     | float, int, float,      | | Create VTK files that contain a slice (cut or average) of the full domain.                     |
|                | string                  | | the "N" of the entry name is an integer that identify each slice, starting from n=1            |
|                |                         | | 1st parameter: Time interval between each slice vtk file                                       |
|                |                         | | 2nd parameter: plane of the slice. 0=(x2,x3) slice, 1=(x1,x3), 2=(x1,x2)                       |
|                |                         | | 3rd parameter: localisation of the slice (when the slice is an average, this parameter only    |
|                |                         | |                affect the localisation of the slice in the produced vtk file                   |
|                |                         | | 4th parameter: slice type. Can be "cut" (for a slice of the full domain) or "average" (for an  |
|                |                         | | average along the direction given by the second parameter). NB: "average" performs a naive     |
|                |                         | | point average, without any consideration on the cell volumes/areas.                            |
|                |                         | | NB2: this feature is in beta, and sometimes fail with some MPI implementations.                |
+----------------+-------------------------+--------------------------------------------------------------------------------------------------+
| xdmf           | float                   | | Time interval between xdmf outputs, in code units (requires Idefix to be configured with HDF5) |
|                |                         | | If negative, periodic xdmf outputs are disabled.                                               |
+----------------+-------------------------+--------------------------------------------------------------------------------------------------+
| xdmf_dir       | string                  | | directory for xdmf file outputs. Default to "./"                                               |
|                |                         | | The directory is automatically created if it does not exist.                                   |
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
|                |                         | | are then written as new variables in vtk and/or xdmf  outputs.                                 |
+----------------+-------------------------+--------------------------------------------------------------------------------------------------+

.. note::
    Even if dumps are not mentionned in your input file (and are therefore disabled), dump files are still produced when *Idefix* captures a signal
    (see :ref:`signalHandling`) or when ``max_runtime`` is set and reached.


.. _dustSection:

``Dust`` section
----------------------

This section describes the dust fields computed using a zero pressure gas approximation (see :ref:`dustModule`).

+----------------+-------------------------+---------------------------------------------------------------------------------------------+
|  Entry name    | Parameter type          | Comment                                                                                     |
+================+=========================+=============================================================================================+
| nSpecies       | integer                 | | Number of dust species to solve                                                           |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| drag           | string, float, ...      | | The first parameter describe the drag type. Possible values are: ``gamma``, ``tau``,      |
|                |                         | | ``size`` and ``userdef``.                                                                 |
|                |                         | | The remaining parameters give the drag parameter :math:`\beta_i` for each dust specie.    |
|                |                         | | (see :ref:`dustModule`). *Idefix* expects to have as many drag parameters as there are    |
|                |                         | | dust species.                                                                             |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| drag_feedback  | bool                    | | (optionnal) whether the gas feedback is enabled (default true).                           |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
