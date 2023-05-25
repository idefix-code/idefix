.. _dustModule:

Dust fluid module
=========================

Equations
---------
The dust module is designed to treat dust grains as a zero-pressure gas coupled to the gas by a linear drag force on the velocity.
*Idefix* can handle as many dust species as one wants, each with different coupling constants. For each dust specie :math:`i`, the code solve the dust evolution equation

.. math::

    \partial_t \rho_i+\mathbf{\nabla}\cdot \rho \mathbf{v}_i&=0

    \partial_t \rho_i \mathbf{v}_i + \nabla\cdot \rho \mathbf{v}_i\mathbf{v}_i&=\gamma_i \rho_i \rho (\mathbf{v}-\mathbf{v}_i)-\rho_i\mathbf{\nabla}\psi_G-\rho_i\mathbf{g}


where :math:`\rho_i` and :math:`\rho` are the dust and gas densities, :math:`\mathbf{v}_i` and :math:`\mathbf{v}` are the dust and gas velocities and :math:`\gamma_i` is the drag coefficient
between the dust and the gas. Note that by construction, all of the source terms (gravity, rotation, shearing box) of the gas are also automatically applied to each dust specie.

When the drag feedback is enabled (default behaviour), the gas is subject to an additional drag force

.. math::

    \mathbf{f}=\sum_i \gamma_i \rho_i \rho (\mathbf{v}_i-\mathbf{v})


which guarantees total momentum conservation.

Drag CFL condition
-------------------
*Idefix* computes the drag terms with a time-explicit scheme. Hence, an addition CFL constraint arrises because of the drag:

.. math::

    dt < \min(\frac{1}{\sum_i\gamma_i(\rho_i+\rho)})

*Idefix* automatically adjusts the CFL to satisfy this inequality, in addition to the usual CFL condition.

Dust parameters
---------------

The dust module can be enabled adding a block `[Dust]` in your input .ini file. The parameters are as follow:

+----------------+-------------------------+---------------------------------------------------------------------------------------------+
|  Entry name    | Parameter type          | Comment                                                                                     |
+================+=========================+=============================================================================================+
| nSpecies       | integer                 | | Number of dust species to solve                                                           |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| drag           | string, float, ...      | | The first parameter describe the drag type. Possible values are: ``gamma``, ``tau``,      |
|                |                         | | ``size`` and ``userdef``.                                                                 |
|                |                         | | The remaining parameters gives the drag parameter :math:`\beta_i` for each dust specie.   |
|                |                         | | (see below). *Idefix* expect to have as many drag parameters as there are dust species.   |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+
| drag_feedback  | bool                    | | (optionnal) whether the gas feedback is enabled (default true).                           |
+----------------+-------------------------+---------------------------------------------------------------------------------------------+

The drag parameter :math:`\beta_i` above sets the functional form of :math:`\gamma_i(\rho, \rho_i, c_s)` depending on the drag type:

``gamma``:
  This sets :math:`\gamma_i=\beta_i`
``tau``:
  This sets :math:`\gamma_i=1/(\rho \beta_i)` so that :math:`\beta_i` is the (constant) stopping time :math:`\tau_i` of the dust grains.
``size``:
  This sets :math:`\gamma_i=1/(c_s \beta_i)` so that :math:`\tau_i=(\rho \gamma_i)^{-1}=\beta_i/(c_s\rho)` where :math:`c_s` is the gas sound speed.
  It designed to reproduce the behaviour of a fixed size dust grain of constant density that follows either Epstein or Stokes drag law.

.. warning::
  The ``userdef`` drag type is not yet implemented.

Using the dust module
---------------------

Several examples are provided in the :file:`test/Dust` directory. Each dust specie is considered in Idefix as a instance of the `Fluid` class, hence
one can apply the technics used for the gas to each dust specie. Because *Idefix* can handle an arbitrarily number of dust species, each specie is stored
in an instance of `Fluid` and stored in a container (:code:`std::vector dust`) in the `DataBlock`. The same is true for the mirror `DataBlockHost`: the
dust primitive variable are all stored in :code:`std::vector dustVc` . For instance, initialising
a single dust specie is done as follow:

.. code-block:: c++


    void Setup::InitFlow(DataBlock &data) {
      // Create a host copy
      DataBlockHost d(data);

      for(int k = 0; k < d.np_tot[KDIR] ; k++) {
          for(int j = 0; j < d.np_tot[JDIR] ; j++) {
              for(int i = 0; i < d.np_tot[IDIR] ; i++) {

                  d.Vc(RHO,k,j,i) = 1.0;            // Set the gas density to 1
                  d.dustVc[0](RHO,k,j,i) = 1.0;     // Set first dust specie density to 1

                  d.Vc(VX1,k,j,i) = 1;              // Set the gas velocity to 1
                  d.dustVc[0](VX1,k,j,i) = 0.0;     // Set the dust velocity to 0

              }
          }
      }

      // Send it all, if needed
      d.SyncToDevice();
    }



All of the dust fields are automatically outputed in the dump and vtk outputs created by *Idefix*.
