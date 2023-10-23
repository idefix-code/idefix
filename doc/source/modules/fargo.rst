.. _fargoModule:

Orbital advection (aka "FARGO") module
=======================================

The Fargo module implements the orbital advection algorithm, to accelerate the computation of flows
subject to a dominant axisymmetric azimuthal velocity component. This is in particular the case
in cold astrophysical discs, for which this method was invented.

The Fargo module implemented in *Idefix* follows the algorithm presented by Mignone et al., A&A 545, 152 (2012).
It can be used in spherical or polar geometries, and in cartesian geometry when combined with the shearing box
boundary conditions.

.. note::
    By default, Fargo uses a piecewise linear advection operator. One can enable
    a piecewise parabolic reconstruction method (ppm) setting ``Idefix_HIGH_ORDER_FARGO``
    to ``ON`` in ``cmake`` configuration. This option is automatically enabled if the
    reconstruction order of the main scheme relies on Limo3 or PPM, for coherence.

The main input of the Fargo module is the fargo velocity, that is, the azimuthal velocity which will be used
as the mean advection velocity in the Fargo scheme. By construction, this velocity is axisymmetric, therefore,
the Fargo velocity field depends only on the coordinates perpendicular to that direction:

* in cartesian geometry, the Fargo velocity depends on :math:`(x,z)` which corresponds to :math:`(x_1,x_3)`
* in cylindrical geometry, the Fargo velocity depends on :math:`(R,z)` which corresponds to :math:`(x_1,x_3)`
* in spherical geometry, the Fargo velocity depends on :math:`(r,\theta)` which corresponds to :math:`(x_1,x_2)`

The Fargo velocity is usually specified enrolling a user-function (see :ref:`functionEnrollment`)
which fills the 2D Fargo velocity arrays (see examples in `test/HD/FargoPlanet` and `test/MHD/FargoMHDSpherical`).
In this case, one should set the `velocity` parameter to `userdef` in the Fargo block of the input file (see :ref:`fargoSection`).
It is also possible to initialse the Fargo velocity field automatically when using the shearing box module in cartesian geometry.
In that case, set the `velocity` parameter to `shearingbox`, and do not forget to enable the `shearingbox` module in the `Hydro` block!

The Fargo module in *Idefix* is fully parallelised. This means that one can have a MPI domain decomposition in any spatial direction, *including*
the azimuthal direction. This last option, however, comes at a cost: the Fargo module needs to exchange a large amount of data with
neighbouring processes for the Fargo advection shift. To limit this communication overhead, Fargo uses a parameter to limit the maximum
number of azimuthal cells over which it will shift the domain at each time step. This optional parameter `maxShift` is by default set to 10.
If it is too small for your setup (i.e. in a case of a very large timestep compared to the mean advection CFL), *Idefix* will stop and tell
you to increase your `maxShift` parameter in the input file. Hence, the user has normally no reason to modify this parameter *a priori*.
