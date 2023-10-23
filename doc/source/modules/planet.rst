.. _planetModule:

Planet module
==============

The planet module allows one to model one or several planets interacting with the gas and other planets. One can force the planet to be on a fixed orbit (in which case
migration and planet-planet interaction are ignored) or on the contrary integrate the equation of motion of the planet(s) simulatenously with the gas. All of the parameters
used by the Planet module are set in the ``[Planet]`` block of the input file and read at runtime.

Common parameters for every planets:

+------------------------+-----------------------+-----------------------------------------------------------------------------------------------------------+
|  Entry name            | Parameter type        | Comment                                                                                                   |
+========================+=======================+===========================================================================================================+
| smoothing              | string, float, float  | | Expression of the planet potential (``plummer`` or ``polynomial``)                                      |
|                        |                       | | and smoothing length for the regularisation of the planet potential around the planet location.         |
|                        |                       | | Mandatory                                                                                               |
+------------------------+-----------------------+-----------------------------------------------------------------------------------------------------------+
| feelDisk               | bool                  | | Whether the planet feels the gravitational forces due to the gas. This parameter enables migration.     |
|                        |                       | | Mandatory                                                                                               |
+------------------------+-----------------------+-----------------------------------------------------------------------------------------------------------+
| feelPlanets            | bool                  | | Whether the planet interacts gravitationnaly with other planets.                                        |
|                        |                       | | Mandatory                                                                                               |
+------------------------+-----------------------+-----------------------------------------------------------------------------------------------------------+
| integrator             | (string)              | | Type of time integrator. Can be ``analytical`` (planets are on fixed orbits) or ``rk4`` (numerical      |
|                        |                       | | integration).                                                                                           |
|                        |                       | | default: ``rk4``                                                                                        |
+------------------------+-----------------------+-----------------------------------------------------------------------------------------------------------+
| indirectPlanets        | (bool)                | | Include indirect planet term arising from the acceleration of the star by the planet.                   |
|                        |                       | | default: ``true``                                                                                       |
+------------------------+-----------------------+-----------------------------------------------------------------------------------------------------------+
| torqueNormalization    | (float)               | | Add a multiplicative factor in the gravitational torque computation                                     |
|                        |                       | | when updating the velocity components of the planet(s) if feelDisk=``true``.                            |
|                        |                       | | default: ``1.0``                                                                                        |
+------------------------+-----------------------+-----------------------------------------------------------------------------------------------------------+
| hillCut                | (bool)                | | Exclude Hill's sphere (with radius :math:`r_\mathrm{h}`) when computing the gravitational torque        |
|                        |                       | | exerted by the gas on the planet. More specifically, we extract the material below                      |
|                        |                       | | :math:`0.5r_\mathrm{h}` and use a :math:`\mathrm{sin}^2` transition until :math:`r_\mathrm{h}`.         |
|                        |                       | | default: ``false``                                                                                      |
+------------------------+-----------------------+-----------------------------------------------------------------------------------------------------------+
| masstaper              | (float)               | | Whether the planet masses are progressively increased or set initially to their values. If set to ``0``,|
|                        |                       | | the masstaper is disabled (default). Otherwise, the value sets the time needed to reach the final       |
|                        |                       | | planets mass.                                                                                           |
+------------------------+-----------------------+-----------------------------------------------------------------------------------------------------------+


Parameters specific to each planet (each entry is followed by a list of values for each planet):

+-----------------------+---------------------+----------------------------------------------------------------------------------------------------------------+
|  Entry name           | Parameter type      | Comment                                                                                                        |
+=======================+=====================+================================================================================================================+
| planetToPrimary       | float, (float, ...) | | Mass ratio between the planet and the central star.                                                          |
|                       |                     | | mandatory                                                                                                    |
+-----------------------+---------------------+----------------------------------------------------------------------------------------------------------------+
| initialDistance       | float, (float, ...) | | Planet initial distance :math:`\displaystyle\frac{x_\mathrm{p}}{1+\mathrm{initialEccentricity}}`.            |
|                       |                     | | mandatory                                                                                                    |
+-----------------------+---------------------+----------------------------------------------------------------------------------------------------------------+
| initialEccentricity   | (float, float, ...) | | Planet initial eccentricity.                                                                                 |
|                       |                     | | default: ``0``                                                                                               |
+-----------------------+---------------------+----------------------------------------------------------------------------------------------------------------+
| initialInclination    | (float, float, ...) | | Planet initial inclination.                                                                                  |
|                       |                     | | default: ``0``                                                                                               |
+-----------------------+---------------------+----------------------------------------------------------------------------------------------------------------+
| tOffset               | (float, float, ...) | | Are the planet(s) directly put in the simulations at t = 0 (default)?                                        |
|                       |                     | | Otherwise, the value sets the time to wait before introducing the planet(s).                                 |
|                       |                     | | default: ``0``                                                                                               |
+-----------------------+---------------------+----------------------------------------------------------------------------------------------------------------+

.. note::
  Note that ``integrator=analytical`` is incompatible with ``feelDisk=true``, ``feelPlanet=true``, ``initialEccentricity!=0``, ``initialInclination!=0``.

``smoothing`` takes 3 values. The first one gives the shape of the planet potential. The last two values (noted here ``a`` and ``b``) give the expression of the smoothing length, of the form :math:`\epsilon=ax_1^b`. In planet-disk interaction modeling, we usually choose ``a = aspect_ratio*thickness_smoothing`` and ``b = flaring_index``, with ``thickness_smoothing = 0.6`` typically in 2D (see, e.g., Masset & Benitez-Llambay, ApJ 817, 19 (2016)).

* A ``plummer`` potential has the form:
.. math::
    \Phi_\mathrm{p}=\displaystyle-\frac{\mathrm{planetToPrimary}}{\sqrt{D_\mathrm{p}^2+\epsilon^2}}
in cartesian coordinates :math:`(x,y,z)`, for a planet located at :math:`(x_\mathrm{p},y_\mathrm{p},z_\mathrm{p})` and with :math:`D_\mathrm{p}=\sqrt{(x-x_\mathrm{p})^2+(y-y_\mathrm{p})^2+(z-z_\mathrm{p})^2}`.

* A ``polynomial`` potential has the form:
.. math::
    \Phi_\mathrm{p}=\left\{
        \begin{array}{ll}
        \displaystyle-\mathrm{planetToPrimary}(\frac{D_\mathrm{p}}{\epsilon}^4 - 2\frac{D_\mathrm{p}}{\epsilon}^3 + 2\frac{D_\mathrm{p}}{\epsilon}) & \mbox{if}~D_\mathrm{p}<\epsilon \\
        \displaystyle-\frac{\mathrm{planetToPrimary}}{D_\mathrm{p}} & \mbox{else}
        \end{array}
        \right.

.. note::
  In 3D, it is possible to simulate only half of the disk by setting one of the boundary in X3 to 0 (in polar coordinates) or X2 to :math:`\pi/2` (in spherical coordinates).
  In this case, we correct the computation of the gravitational torque exerted by the gas onto the planet(s) by removing the contribution of the vertical component of the torque and doubling
  the horizontal contributions (hence assuming midplane symmetry).
