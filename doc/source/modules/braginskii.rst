.. _braginskiiModule:

Braginskii module
===================

Equations solved and method
---------------------------

The ``Braginskii`` module implements the anisotropic heat and momentum fluxes specific
to weakly collisional, magnetised plasma like the intracluster medium
for which this module was originally designed.
Physically, it computes the modified Fourier law :math:`q_\mathrm{B}` and
the modified viscous stress tensor :math:`\Pi_\mathrm{B}`:

:math:`q_\mathrm{B} = \kappa_\parallel \left(\hat{b}\cdot \nabla T\right) \hat{b}`
:math:`\Pi_\mathrm{B} = - \left( p_\perp - p_\parallel \right)  \left( \hat{b} \hat{b} - \frac{1}{3} \vec{I} \right)`

(where :math:`\kappa_\parallel` is the parallel thermal conductivity in code units.
The pressure anisotropy can be computed from the following closure:
:math:`p_\perp - p_\parallel = 3\mu_\mathrm{B} \left(\hat{b}\hat{b}:\nabla v - \frac{1}{3} \nabla\cdot v \right)`.

The anisotropic heat flux from the ``Braginskii`` module is implemented in *Idefix*
according to the centered asymmetric scheme described in Sharma & Hammet (2007, Section 2.1).
The Braginskii viscous stress tensor is implemented with the same scheme,
though adapated to vector quantities.

.. note::
    By default, the Braginskii module
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
|                |                         | | the solver  name (e.g. ``PCG`` or ``PBIGCSTAB`` ).                                        |
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

The Braginskii module in *Idefix* is fully parallelised. This means that one can have a MPI domain decomposition in any spatial direction
either on CPU or GPU.
