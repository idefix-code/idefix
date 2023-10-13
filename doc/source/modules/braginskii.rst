.. _braginskiiModule:

Braginskii module
===================

Equations solved and methods
---------------------------

The ``Braginskii`` module implements the anisotropic heat and momentum fluxes specific
to weakly collisional, magnetised plasma like the intracluster medium
for which this module was originally designed.
In practice, it computes the modified Fourier law :math:`q_\mathrm{B}` and
the modified viscous stress tensor :math:`\Pi_\mathrm{B}` (Latter & Kunz 2012):

:math:`q_\mathrm{B} = -\kappa_\parallel \left(\hat{b}\cdot \nabla T\right) \hat{b}`

:math:`\Pi_\mathrm{B} = - \left( p_\perp - p_\parallel \right)  \left( \hat{b} \hat{b} - \frac{1}{3} I \right)`

where
:math:`T` and :math:`v` are respectively the temperature and the velocity of the fluid,
:math:`p_\parallel` and :math:`p_\perp` the pressure in the direction parallel and
perpendicular to the local magnetic field.
:math:`\kappa_\parallel` and :math:`\mu_\parallel` are the parallel thermal conductivity
and dynamic viscosity respectively, while
:math:`\hat{b}` is the unit vector in the direction of the local magnetic field
and :math:`I` the identity.
The pressure anisotropy can be computed from the following closure:
:math:`p_\perp - p_\parallel = 3\mu_\mathrm{B} \left(\hat{b}\hat{b}:\nabla v - \frac{1}{3} \nabla\cdot v \right)` (Schekochihin et al. 2010).

The anisotropic heat flux from the ``Braginskii`` module is implemented in *Idefix*
according to the centered asymmetric scheme described in Sharma & Hammett (2007, Section 2.1).
The Braginskii viscous stress tensor is implemented with the same scheme,
though adapted to vector quantities.

.. note::
    By default, the Braginskii module computes the transverse heat flux terms at the right
    cell interface by a simple arithmetic average
    (Eq. (6)-(7) from Sharma & Hammett 2007).
    However in the same paper, the authors showed that this implementation can lead to
    unphysical heat flux from high to low temperature regions.
    So we implemented slope limiters for the computation of these transverse heat fluxes,
    as described in Eq. (17) from Sharma & Hammett (2007).
    Only the van Leer and the Monotonized Central (MC) limiters are available
    since the minmod limiter has been found to be too diffusive.
    Transverse viscous fluxes are limited in the same manner,
    as described by ZuHone et al. (2015) in their fifth footnote.

.. note::
    Perpendicular heat flux
    :math:`q_\perp = -\kappa_\perp \left[ \nabla T - \left(\hat{b}\cdot \nabla T\right) \hat{b}\right]`
    can also be taken into account
    (see section :ref:`braginskiiParameterSection`),
    but perpendicular viscous flux cannot in the current implementation.
    However, both Braginskii operators can be used *in addition* to the classical
    Fourier law and viscosity. Modeling both isotropic and anisotropic diffusion
    should be equivalent to having different
    parallel and perpendicular diffusivities.

The main output of the ``Braginskii`` module is the addition of the anisotropic heat and/or
momentum flux terms to the other various fluxes.
The computational cost of the ``Braginskii`` module is therefore one
of a parabolic term, meaning that
the associated Courant-Friedrichs-Lewy (CFL)
condition can be stiff in case of very diffusive plasma.
However, like other parabolic operators in *Idefix*, a Runge-Kutta-Legendre super time-stepping
scheme is available to speed-up the computation (with respect to a fully explicit integration)
of the Braginskii heat flux and viscosity.

..
  Please refer to Section 2.8 from Lesur et al.
  for more details on the this super time-stepping scheme.

.. _braginskiiParameterSection:

Main parameters of the module
-----------------------------

The ``Braginskii`` module can be enabled adding one or two lines in the ``[Hydro]`` section
starting with the keyword
`bragTDiffusion` or/and *bragViscosity*. The following table summarises the different options
associated to the activation of the Braginskii module:

+--------+-----------------------+-------------------------+---------------------------------------------------------------------------------------+
| Column |  Entry name           | Parameter type          | Comment                                                                               |
+========+=======================+=========================+=======================================================================================+
| 0      |  bragModule           | string                  | | Activates Braginskii diffusion. Can be ``bragTDiffusion`` or ``bragViscosity``.     |
+--------+-----------------------+-------------------------+---------------------------------------------------------------------------------------+
| 1      | integration           | string                  | | Specifies the type of scheme to be used to integrate the parabolic term.            |
|        |                       |                         | | Can be ``rkl`` or ``explicit``.                                                     |
+--------+-----------------------+-------------------------+---------------------------------------------------------------------------------------+
| 2      | slope limiter         | string                  | | Choose the type of limiter to be used to compute anisotropic transverse flux terms. |
|        |                       |                         | | Can be ``mc``, ``vanleer`` or ``nolimiter``.                                        |
+--------+-----------------------+-------------------------+---------------------------------------------------------------------------------------+
| 3      | diffusivity type      | string                  | | Specifies the type of diffusivity wanted. Can be ``constant`` or ``userdef``.       |
+--------+-----------------------+-------------------------+---------------------------------------------------------------------------------------+
| 4      | parallel diffusivity  | real                    | | Mandatory if the diffusivity type is ``constant``. Not needed otherwise.            |
|        |                       |                         | | Value of the parallel diffusivity. Should be a real number.                         |
+--------+-----------------------+-------------------------+---------------------------------------------------------------------------------------+
| 5      | normal diffusivity    | real                    | | When bragModule ``bragTDiffusion`` and diffusivity type ``constant``,               |
|        |                       |                         | | value of the normal diffusivity. Should be a real number.                           |
+--------+-----------------------+-------------------------+---------------------------------------------------------------------------------------+

Numerical checks
---------------

In Cartesian geometry, the ``Braginskii`` module has been tested with many setups
and in all configurations of magnetic polarisation:
growth rates of local instabilities (see the MTI test inspired from Parrish et al. 2012)
and damped rates of eigenmodes of the corresponding Braginskii operator (tests not included).
In Cylindrical/Polar geometry, only the anisotropic heat conduction has been tested
with numerical measurements of the damped rates of its eigenmodes, in all directions
(tests not included).
In Spherical geometry, both Braginskii operators have been validated by measuring the damped rates
of their eigenmodes for a purely radial and purely azimuthal magnetic polarisation
(tests included except the viscosity with an azimuthal magnetic field).

In conclusion, the ``Braginskii`` operators have been thoroughly tested in Cartesian geometry.
The same goes for the anisotropic heat flux in Cylindrical/Polar geometry while
the anisotropic viscosity has *never* been tested in this geometry.
In spherical geometry, both ``Braginskii`` operators have been partially validated
(diffusion along the polar axis has not been directly tested).
