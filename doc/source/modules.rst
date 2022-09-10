=====================
*Idefix* modules
=====================

In this section, you will find a more detailed documentation about each module that can be used in
*Idefix*

:ref:`fargoModule`
  The orbital advection module, speeds up the computation of flows dominated by an azimuthal motion (such as discs).

:ref:`selfGravityModule`
  The self-gravity computation module, handles the impact of the gas distribution on its own dynamic when massive
  enough (as in a core collapse).

:ref:`gridCoarseningModule`
  The grid coarsening module, that allows to derefine the grid in particular locations to speed up the computation.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules/fargo.rst
   modules/selfGravity.rst
   modules/gridCoarsening.rst
