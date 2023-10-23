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

When the gas follows an ideal or user-defined equation of state, the dust-gas drag leads to friction heating. Because the dust component does not posess an internal energy, *Idefix*
assumes that all of the friction heating is deposited in the gas internal energy as an addition heating term :math:`Q_d`:

.. math::

    Q_d = \sum_i \gamma_i \rho_i \rho (\mathbf{v}-\mathbf{v}_i)\cdot\mathbf{v}_i

When the drag feedback is enabled (default behaviour), the gas is subject to an additional drag force

.. math::

    \mathbf{f}=\sum_i \gamma_i \rho_i \rho (\mathbf{v}_i-\mathbf{v})

which guarantees total momentum conservation.

.. note::

    The heating associated to the feedback does not require any additional source term in the energy equation as *Idefix* conserves the total gas energy by design.



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
  This sets :math:`\gamma_i=c_s/\beta_i` so that :math:`\tau_i=(\rho \gamma_i)^{-1}=\beta_i/(\rho c_s)` where :math:`c_s` is the gas sound speed.
  It designed to reproduce the behaviour of a fixed size dust grain of constant density that follows either Epstein or Stokes drag law:

    Epstein drag:
      In the Epstein regime, :math:`\beta_i=\rho_s a` where :math:`\rho_s` is the solid density and :math:`a` is the solid size.
    Stokes drag:
      In the Stokes regime, :math:`\beta_i=4\rho_s a^2/(9\lambda_\mathrm{mfp})` where :math:`\lambda_\mathrm{mfp}` is the gas mean free path.

``userdef``:
  This allows the user to define a specialized drag law, that is a function :math:`\gamma_i(\rho, \rho_i, c_s, \beta_i)`. In this case, a user-defined
  function should be enrolled by each drag instance (see the example in `test/Dust/FargoPlanet`). Note that the entry ``drag`` in your
  input file should still contain a list of :math:`\beta_i`.


.. warning::
  The drag force assumes that the gas density field ``hydro->Vc(RHO...)`` is a volumic density. For 2D problems assuming
  a razor-thin geometry, this assumption is incorrect since ``hydro->Vc(RHO...)`` is a surface density. In this case,
  the user has to define a specific drag law since *Idefix* has no way to guess how to convert the surface density to
  the volumic density (see example in `test/Dust/FargoPlanet`).



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
