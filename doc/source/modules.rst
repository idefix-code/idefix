=====================
*Idefix* modules
=====================

In this section, you will find a more detailed documentation about each module that can be used in
*Idefix*

:ref:`fargoModule`
  The orbital advection module, speeds up the computation of flows dominated by an azimuthal motion (such as discs).

:ref:`planetModule`
  The planet module, which treats the planet-disk interaction and planet-planet interaction.

:ref:`dustModule`
  The dust module, modeling dust grains as a zero-pressure gas.

:ref:`eosModule`
  The custom equation of state module, allowing the user to define its own equation of state.

:ref:`selfGravityModule`
  The self-gravity computation module, handles the impact of the gas distribution on its own dynamic when massive
  enough (as in a core collapse).

:ref:`braginskiiModule`
  The Braginskii module, models the anisotropic flux of heat and momentum
  taking place in weakly collisional, magnetised plasma (like the intracluster medium).

:ref:`gridCoarseningModule`
  The grid coarsening module, that allows to derefine the grid in particular locations to speed up the computation.

:ref:`pydefixModule`
  The Python interaction module, that allows Idefix to interact directly with a python interpreter.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules/fargo.rst
   modules/planet.rst
   modules/dust.rst
   modules/eos.rst
   modules/selfGravity.rst
   modules/braginskii.rst
   modules/gridCoarsening.rst
   modules/pydefix.rst
