.. _selfGravityModule:

Self Gravity module
===================

Equations solved and method
---------------------------

The ``SelfGravity`` module implements the computation of a self-gravitational potential generated
by the gas distribution. Physically, it solves for the potential :math:`\psi_{SG}` that satisfies the
Poisson equation

:math:`\Delta \psi_{SG}=4\pi G_c\,\rho` (where :math:`G_c` is the gravitational constant in code
units)

This computation becomes relevant when the said gas distribution
is massive enough to influence its own dynamic. This is the case for example in young protoplanetary
discs, which are subject to gravitational instabilities and for which this module was originally designed.

The ``SelfGravity`` module implemented in *Idefix* follows the algorithm of the BIConjugate Gradient
STABilized (BICGSTAB) method, a common iterative method designed to solve sparse linear systems (see for example
Saad, Y. (2003). Iterative methods for sparse linear systems. Society for Industrial and Applied Mathematics.).
This specific method allows to solve the Poisson equation in any dimensions and any geometries.

.. note::
    By default, the selfGravity module uses the BICGSTAB method to invert the Poisson equation.
    However, an optional preconditionning feature (PBICGSTAB) is
    implemented to handle particularly irregular and inhomogeneous grids, that may prevent the
    convergence of the classic BICGSTAB. Note also that a simpler (yet very slow) Jacobi method
    has been left for debug purpose. The user can also try the conjugate gradient and minimal residual
    methods which have been tested successfully and are faster than BICGSTAB for some problems/grids.

The main output of the ``SelfGravity`` module is the addition of the self-gravitational potential inferred from the
gas distribution to the various sources of gravitational potential. At the beginning of every (M)HD step, the module is called to compute
the potential due to the mass distribution at the given time. The potential computed by the ``SelfGravity`` module
is then added to the other potential sources by the ``Gravity`` class (planets, central and user-defined potential).
The computational cost of one selfGravity computation depends on the complexity of the
gravitational field and the speed at which the density field evolves.

Main parameters of the module
-----------------------------

The ``SelfGravity`` module is a submodule of the ``Gravity`` module to compute the gravitationnal potential. Hence, it is enabled
by adding ``selfgravity`` to the ``potential`` in the ``[Gravity]`` block (see :ref:`gravitySection`). The parameters specific to self gravity are to be
set in a dedicated ``[SelfGravity]`` block:

+----------------+-------------------------+---------------------------------------------------------------------------------------------+
|  Entry name    | Parameter type          | Comment                                                                                     |
+================+=========================+=============================================================================================+
| solver         | string                  | | Specifies which solver should be used. Can be ``Jacobi``, ``CG``, ``MINRES``, ``BICGSTAB``|
|                |                         | | which corresponds to Jacobin, conjugate gradient, Minimal residual or bi-conjugate        |
|                |                         | | stabilised method. Note that a preconditionned version is available adding a ``P`` to     |
|                |                         | | the solver name (e.g. ``PCG`` or ``PBICGSTAB``).                                         |
|                |                         | | The FFT solver (``FFT``) is also available for Cartesian geometries with periodic          |
|                |                         | | boundaries (see :ref:`fftSolverSection`).                                                |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| targetError    | real                    | | Set the error allowed in the residual :math:`r=\Delta\psi_{SG}/(4\pi G_c)-\rho`. The error|
|                |                         | | computation is based on a L2 norm. Default is 1e-2.                                       |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| maxIter        | int                     | | Set the maximum number of iterations allowed to the solver to reach convergence. Default  |
|                |                         | | is 1000.                                                                                  |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| skip           | int                     | | Set the number of integration cycles between each computation of self-gravity potential.  |
|                |                         | | Default is 1 (i.e. self-gravity is computed at every cycle).                              |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+




Boundary conditions on self-gravitating potential
--------------------------------------------------

Solving for the self-gravitating potential require the definition of a boundary condition for said potential. This is done with the entries
``boundary-Xn-dir`` of the ``[SelfGravity]`` block, where ``n`` can be 1, 2 or 3 and is the direction for the boundary condition and ``dir`` can be ``beg`` or ``end`` and
indicates the side of the boundary.

The boundary conditions can be following

+-----------------------+------------------------------------------------------------------------------------------------------------------+
| ``boundary-Xn-dir``   | Comment                                                                                                          |
+=======================+==================================================================================================================+
| nullpot               | Zero potential boundary conditions. The potential is set to zero in the ghost cells.                             |
+-----------------------+------------------------------------------------------------------------------------------------------------------+
| nullgrad              | Zero gradient boundary conditions. The potential in the last active cell is copied in the ghost cells.           |
+-----------------------+------------------------------------------------------------------------------------------------------------------+
| periodic              | Periodic boundary conditions. The potential is copied between beg and end sides of the boundary.                 |
+-----------------------+------------------------------------------------------------------------------------------------------------------+
| axis                  | | Axis boundary condition. Should be used in spherical coordinate in the X2 direction when the domain starts/stop|
|                       | | on the axis.                                                                                                   |
+-----------------------+------------------------------------------------------------------------------------------------------------------+
| origin                | | In spherical coordinates, artificially extends the grid used to compute the potential close to R=0.            |
|                       | | should only be used in X1-beg direction.                                                                       |
+-----------------------+------------------------------------------------------------------------------------------------------------------+
| userdef               | |User-defined boundary conditions. The boundary condition function should be enrolled in the setup constructor   |
|                       | | (see :ref:`userdefBoundaries`).                                                                                |
+-----------------------+------------------------------------------------------------------------------------------------------------------+

.. note::
    The method in fully periodic setups requires the removal of the mean gas density
    before solving Poisson equation. This is done automatically if all of the self-gravity boundaries are set to ``periodic``.
    Hence, make sure to specify all self-gravity boundary conditions as periodic for such setups, otherwise the solver will
    fail to converge.


.. _fftSolverSection:

FFT Solver
----------

The FFT-based solver provides an alternative to iterative methods for solving the Poisson equation.
It computes the gravitational potential in Fourier space, where Poisson's equation becomes a simple
algebraic relation: :math:`\hat{\psi} = -4\pi G_c \hat{\rho} / k^2` (for :math:`k \neq 0`).
This approach is significantly faster than iterative methods for uniform Cartesian grids.

**Requirements and limitations:**

- **Geometry**: Cartesian only (not compatible with cylindrical or spherical geometries)
- **Dimensionality**: Only 3D is supported, as the FFT solver relies on 3D FFTs for efficiency
- **Boundary conditions**: Must be periodic in all active dimensions
- **Grid constraints**:
  - Even number of grid points required in the X1 direction
  - No domain decomposition allowed in the X1 and X2 directions
  - The number of grid points in X2 and X3 must be divisible by the number of MPI processes (due to global transpositions)
- **Compilation**: Requires enabling FFT support with ``-DIdefix_FFT=ON`` during CMake configuration

**Parallelization:**

The FFT solver uses a slab decomposition strategy for MPI parallelization, where the domain can be
decomposed along the X2 and X3 directions. The solver leverages the KokkosFFT library for efficient
FFT computations on both CPU and GPU architectures.
