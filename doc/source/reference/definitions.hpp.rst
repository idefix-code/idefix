=======================================
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
