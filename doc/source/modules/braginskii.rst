Braginskii module
===================

Equations solved and methods
---------------------------

The ``Braginskii`` module implements the anisotropic heat and momentum fluxes specific
to weakly collisional, magnetised plasma like the intracluster medium
for which this module was originally designed.
Physically, it computes the modified Fourier law :math:`q_\mathrm{B}` and
the modified viscous stress tensor :math:`\Pi_\mathrm{B}` (Parrish et al. 2012):

:math:`q_\mathrm{B} = -\kappa_\parallel \left(\hat{b}\cdot \nabla T\right) \hat{b}`

:math:`\Pi_\mathrm{B} = - \left( p_\perp - p_\parallel \right)  \left( \hat{b} \hat{b} - \frac{1}{3} I \right)`

where
:math:`T` and :math:`v` are respectively the temperature and the velocity of the fluid,
:math:`p_\parallel` and :math:`p_\perp` the pressure in the direction parallel and
perpendicular to the local magnetic field.
:math:`\kappa_\parallel` and :math:`\mu_\parallel` are the parallel thermal conductivity
and dynamic viscosity respectively,
:math:`\hat{b}` the unit vector in the direction of the local magnetic field
and :math:`I` the identity.
The pressure anisotropy can be computed from the following closure:
:math:`p_\perp - p_\parallel = 3\mu_\mathrm{B} \left(\hat{b}\hat{b}:\nabla v - \frac{1}{3} \nabla\cdot v \right)`.

The anisotropic heat flux from the ``Braginskii`` module is implemented in *Idefix*
according to the centered asymmetric scheme described in Sharma & Hammett (2007, Section 2.1).
The Braginskii viscous stress tensor is implemented with the same scheme,
though adapted to vector quantities.

By default, the Braginskii module computes the transverse heat flux terms at the right
cell interface by a simple arithmetic average
(Eq. (6)-(7) from Sharma & Hammett 2007).
However in the same paper, the authors showed that this implementation can lead to
unphysical heat flux from high to low temperature regions.
We therefore implemented slope limiters for the computation of these transverse fluxes,
as described in Eq. (17) from Sharma & Hammett (2007).
Only the van Leer and the Monotonized Central (MC) limiters are available
since minmod has been found to be too diffusive.
Transverse viscous fluxes are limited in the same manner,
as described by ZuHone et al. (2015) in their fifth footnote.

..
   .. note::
    this::By default, the Braginskii module computes the transverse heat flux terms at the right
    cell interface by a simple arithmetic average
    (Eq. (6)-(7) from Sharma & Hammett 2007).
    However in the same paper, the authors showed that this implementation can lead to
    unphysical heat flux from high to low temperature regions.
    We therefore implemented slope limiters for the computation of these transverse fluxes,
    as described in Eq. (17) from Sharma & Hammett (2007).
    Only the van Leer and the Monotonized Central (MC) limiters are available
    since minmod has been found to be too diffusive.
    Transverse viscous fluxes are limited in the same manner,
    as described by ZuHone et al. (2015) in their fifth footnote.

.. note::
    Perpendicular heat flux
    :math:`q_\perp = -\kappa_\perp \left[ \nabla T - \left(\hat{b}\cdot \nabla T\right) \hat{b}\right]`
    can also be taken into account
    (see section :ref:`braginskiiParameterSection`),
    but perpendicular viscous flux cannot in the current implementation.
    However, both Braginskii operators can be used *in addition* to the classical
    Fourier law and viscosity, which should be equivalent to having different
    parallel and perpendicular diffusivities.

The main output of the ``Braginskii`` module is the addition of the anisotropic heat and/or 
momentum flux terms to the other various fluxes.
The computational cost of the ``Braginskii`` module is one of a parabolic term, meaning that
the associated Courant-Friedrichs-Lewy condition can be stiff in case of very diffusive plasma.
However, as other parabolic operators in *Idefix*, a Runge-Kutta-Legendre scheme is available
to speed-up the computation (with respect to an explicit integration)
of the Braginskii heat flux and viscosity.
Please refer to Section 2.8 from Lesur et al.
for more details on the this super time-stepping scheme.

.. _braginskiiParameterSection:
Main parameters of the module
-----------------------------

The ``Braginskii`` module can be enabled adding one or two lines in the ``[Hydro]`` section
starting with the keyword
`bragTDiffusion` or/and *bragViscosity*. The following table summarises the different options
associated to the activation of the Braginskii module:

+----------------+-------------------------+---------------------------------------------------------------------------------------------+
|  Entry name    | Parameter type          | Comment                                                                                     |
+================+=========================+=============================================================================================+
| bragModule     | string                  | | Activates Braginskii diffusion. Can be ``bragTDiffusion`` or ``bragViscosity``            |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| integration sch| string                  | | Specifies the type of scheme to be used to integrate the parabolic term. Can be ``rkl`` or|
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| integration sch| string                  | | Specifies the type of scheme to be used to integrate the parabolic term. Can be ``rkl`` or|
+----------------+-------------------------+---------------------------------------------------------------------------------------------+

..
  +----------------+-------------------------+---------------------------------------------------------------------------------------------+
  Column number |  Entry name    | Parameter type          | Comment                                                                                     |
  +================+=========================+=============================================================================================+
  0 | bragTDiffusion         | string                  | | Activate Braginskii heat diffusion. |
  +----------------+-------------------------+---------------------------------------------------------------------------------------------+
  1 | explicit/rkl           | string                  | | 
  +----------------+-------------------------+---------------------------------------------------------------------------------------------+
  | maxIter        | int                     | | Set the maximum number of iterations allowed to the solver to reach convergence. Default  |
  |                |                         | | is 1000.                                                                                  |
  +----------------+-------------------------+---------------------------------------------------------------------------------------------+
  | skip           | int                     | | Set the number of integration cycles between each computation of self-gravity potential.  |
  |                |                         | | Default is 1 (i.e. self-gravity is computed at every cycle).                              |
  +----------------+-------------------------+---------------------------------------------------------------------------------------------+


..
  this::The Braginskii module in *Idefix* is fully parallelised. This means that one can have a MPI domain decomposition in any spatial direction either on CPU or GPU.
