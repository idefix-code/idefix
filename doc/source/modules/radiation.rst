.. _radiationModule:

Radiation module
=========================

Equations
---------
The radiation module is designed to treat the coupled evolution of a (magneto-)hydrodynamic fluid and a radiative field.
The radiative field is effectively treated as a radiative fluid by solving for the first two moments of the grey radiative transfer equation, the so-called M1 method.
The full set of HD-radiation equation can be written as

.. math::
    \partial_t \rho + \nabla \cdot \left[ \rho \mathbf{v} \right] &= 0,\\
    \partial_t (\rho \mathbf{v}) + \nabla \cdot \left[ \rho \mathbf{v} \otimes \mathbf{v} + \mathcal{P} \mathbb{I} \right] &= -\rho \nabla \psi + \mathbf{G}, \\
    \partial_t E + \nabla \cdot \left[ (E + \mathcal{P}) \mathbf{v} \right] &= -\rho \mathbf{v} \cdot \nabla \psi + c G^0 - \nabla \cdot \mathbf{F}_\mathrm{irr},\\
    \partial_t E_r + \nabla \cdot \mathbf{F}_r &= -\hat{c}\:G^0, \\
    \partial_t \mathbf{F}_r + \nabla \cdot \mathbb{P}_r &= - \hat{c}\:\mathbf{G},

where :math:`E_r`, :math:`\mathbf{F}_r` and :math:`\mathbb{P}_r` are the radiative energy, the radiative energy flux and the radiative energy tensor defined respectively from the direction and frequency-dependent specific intensity, :math:`I_v(t,\mathbf{x},\mathbf{n})` as

.. math::
    E_r &= \frac{1}{c} \int_0^\infty d\nu \oint d\Omega \, I_\nu(t, \mathbf{x}, \mathbf{n}) \\
    F_r^i &= \frac{1}{c} \int_0^\infty d\nu \oint d\Omega \, I_\nu(t, \mathbf{x}, \mathbf{n})\, n^i \\
    P_r^{ij} &= \frac{1}{c} \int_0^\infty d\nu \oint d\Omega \, I_\nu(t, \mathbf{x}, \mathbf{n})\, n^i n^j.


To close the set of the radiation moment equations, we use the closure from Levermore, C. D. 1984, J. Quant. Spectr. Rad. Transf., 31, 149 :

.. math::
    P^{ij}_r = D^{ij} E_r,
with

.. math::
    D^{ij} = \frac{1 - \xi}{2} \delta^{ij} + \frac{3\xi - 1}{2} n^i n^j,
and

.. math::
    \xi = \frac{3 + 4f^2}{5 + 2\sqrt{4 - 3f^2}},

where :math:`\delta^{ij}` is the Kronecker delta symbol, :math:`\mathbf{n}=\mathbf{F_r}/\lVert\mathbf{F_r}\rVert` is a unit vector pointing in the direction of the radiative flux, and :math:`f` is the so-called reduced flux defined as

.. math::
    f\equiv \frac{\lVert\mathbf{F_r}\rVert}{E_r}\le 1.


The source terms in the radiation equations are defined as

.. math::
    G^0 &= \kappa_\mathrm{P} \, \rho \left( E_r - a_R T^4 \right)  \\
    \mathbf{G} &= \chi \, \rho \, (\mathbf{F}_r-4 E_r\mathbf{\beta}).



where :math:`\kappa_\mathrm{P}` is the Planck absorption opacity, :math:`\chi\equiv\kappa_\mathrm{R}+\sigma` where :math:`\kappa_\mathrm{R}` and :math:`\sigma` are the Rosseland absorption and scattering opacities respectively.
The present formulation of the source terms is only valid for non-relativistic sources with :math:`\beta\equiv\frac{||\mathbf{v}||}{c} \ll 1` so that we only take into account a correction of order :math:`\tau \beta` to model the dynamic diffusion limit.
For midly relativistic sources, higher order corrections that are included here are likely to be important.

In addition to the matter-radiation interaction source terms, the user can add an additional heating term :math:`- \nabla \cdot \mathbf{F}_\mathrm{irr}`.
This irradiation source term is not part of the radiation scheme.

.. note::

    The radiation module is incompatible with an isothermal equation of state.



Time Integration
----------------

Radiation is solved using an explicit-implicit method. The hyperbolic part of the radiative equations are solved using the finite-volume Godunov high order schemes of *Idefix*.
The source terms are solved using an implicit method.

Explicit part: Reduced speed of light approximation
+++++++++++++++++++++++++++++++++++++++++++++++++++

In order to reduce the CFL constraint that would arise when dealing with explicit radiation transport at the speed of light, we use the reduced speed of light approximation :math:`\hat{c}`.
The value of  :math:`\hat{c}` is set by the user and should be chosen so that radiation sets the shortest time-scale.

Two explicit solvers are implemented, ``hll_rad`` and ``hllc_rad``, although we exclusively use hll_rad for our applications so that ``hllc_rad`` is not tested for complex setups.

Implicit part:
+++++++++++++++

The source terms are included using an implicit scheme. All of the tests have been carried with the full_implicit solver but two iterative other solvers (fixed_point_rad and fixed_point_gas) described in ??? are also implemented although less tested.

When using the ``full_implicit`` solver the radiation energy and the radiation fluxes are evolved by inverting the following system

.. math::
    M =
    \begin{pmatrix}
    1+\gamma_\mathrm{abs,red} & -4\gamma_\mathrm{abs,red} a_R {T^n}^3 & 0 & 0 & 0 \\
    -\gamma_\mathrm{abs} & C_v+4\gamma_\mathrm{abs} a_R {T^n}^3 & 0 & 0 & 0 \\
    -2\beta_1\gamma_\mathrm{scatt,red}/3 & 0 & 1+\gamma_\mathrm{scatt,red} & 0 & 0 \\
    -2\beta_1\gamma_\mathrm{scatt,red}/3 & 0 & 0 & 1+\gamma_\mathrm{scat,red} & 0 \\
    -2\beta_1\gamma_\mathrm{scatt,red}/3 & 0 & 0 & 0 & 1+\gamma_\mathrm{scatt,red}
    \end{pmatrix}

with the following right hand side

.. math::
    S =
    \begin{pmatrix}
    1+\gamma_1 \\
    -\gamma_2  \\
    -C_1\gamma_3 \\
    -C_2\gamma_3  \\
    -C_3\gamma_3
    \end{pmatrix}

Then, the gas energy and gas momentum are updated to conserve to total energy and total momentum by doing

.. math::


Radiation parameters
---------------

The radiation module can be enabled adding a block `[Rad]` in your input .ini file. The parameters are as follow:

+----------------+-------------------------+----------------------------------------------------------------------------------------------------------------+
|  Entry name    | Parameter type          | Comment                                                                                                        |
+================+=========================+================================================================================================================+
| nFrequencies   | integer                 | | Number of frequency group for radiation.                                                                     |
|                |                         | | The current implementation of *Idefix* allows for only one frquency group but we plan                        |
|                |                         | | to extend our method to several multigroup in the near future                    .                           |
+----------------+-------------------------+----------------------------------------------------------------------------------------------------------------+
| solver_rad     | string                  | | Riemann solver for the hyperbolic part of the radiation equation.                                            |
|                |                         | | Possible values are: ``lfr_rad`` for a radiative Lax–Friedrichs–Rusanov solver                               |
|                |                         | | ``hll_rad`` for the HLL radiative solver introduced in Gonz\'{a}lez et al. (2007)                            |
|                |                         | | ``hllc_rad`` for the HLLC radiative solver introduced in Melon Fuksman \& Mignone (2019)                     |
+----------------+-------------------------+----------------------------------------------------------------------------------------------------------------+
| reduced_c      | float                   | | Ratio of the reduced speed of light to the speed of light.                                                   |
+----------------+-------------------------+----------------------------------------------------------------------------------------------------------------+
| kappa          | string, ...             | | Sets the Planck and Rosseland absorption opacities.                                                          |
|                |                         | | The first parameter sets the absorption opacity type.                                                        |
|                |                         | | For ``constant``, the second and third parameters set a constant Planck and Rosseland absorption opacity     |
|                |                         | | respectively.                                                                                                |
|                |                         | | For ``kramers``, opacity has the form :math:`\kappa=\kappa_0(\frac{\rho}{\rho_0})(\frac{T}{T_0})^{-3.5}`.    |
|                |                         | | The second, third, fourth and fifth parameters set :math:`\kappa_{P,0}`, :math:`\kappa_{R,0}`,               |
|                |                         | | and  :math:`\rho_0` and :math:`T_0` respectively.                                                            |
|                |                         | | For ``usertable``, the absorption opacities are interpolated from a usertable.                               |
|                |                         | | The second parameter sets the dimension of the opacity tables (only 1 or 2 dimensions).                      |
|                |                         | | The third and fourth parameters should be the name of the files containing the Planck and Rosseland          |
|                |                         | | absorption opacity tables respectively.                                                                      |
|                |                         | | For ``userfunc``, absorption opacities are defined from a user defined function. In this case, ``RadSource`` |
|                |                         | | class expects a function to be enrolled with ``Radiation[nFrequencies]::EnrollKappa(KappaFunc)``.            |
+----------------+-------------------------+----------------------------------------------------------------------------------------------------------------+
| xi             | string, ...             | | Sets the scattering opacity.                                                                                 |
|                |                         | | The first parameter sets the scattering opacity type.                                                        |
|                |                         | | For ``constant``, the second parameter sets the constant scattering opacity.                                 |
|                |                         | | For ``usertable``, the scattering opacity is interpolated from a usertable.                                  |
|                |                         | | The second parameter sets the dimension of the opacity table (only 1 or 2 dimensions).                       |
|                |                         | | The third parameter should be the name of the file containing the scattering opacity.                        |
|                |                         | | For ``userfunc``, scattering opacity is defined from a user defined function. In this case, ``RadSource``    |
|                |                         | | class expects a function to be enrolled with ``Radiation[nFrequencies]::EnrollXi(XiFunc)``                   |
+----------------+-------------------------+----------------------------------------------------------------------------------------------------------------+
| source         | string                  | | Implicit method for the matter/radiation interaction source terms.                                           |
|                |                         | | For ``full_implicit`` the source terms are written in the form of a matrix that is inverted analtyically by  |
|                |                         | | linearizing :math:`({T^{n+1})}^4 \approx 4T^{n+1}({T^n})^3-3({T^{n}})^4`. It is the fastest and most robust  |
|                |                         | | method. It is also the most tested.                                                                          |
|                |                         | | For ``fixed_point_rad``, the source terms are solved using a fixed point iterative method where iteration    |
|                |                         | | is made on the computation of the radiative energy. This method is to be used when radiation pressure        |
|                |                         | | dominates over gas pressure.                                                                                 |
|                |                         | | For ``fixed_point_gas``, the source terms are solved using a fixed point iterative method where iteration    |
|                |                         | | is made on the computation of the gas energy. This method is to be used when gas pressure                    |
|                |                         | | dominates over radiation pressure.                                                                           |
+----------------+-------------------------+----------------------------------------------------------------------------------------------------------------+
| irr            | string, ...             | | (optional, default is no irradiation source) Sets the irradiation source.                                    |
|                |                         | | The first parameter sets the type of irradiation source.                                                     |
|                |                         | | Possible values are: ``constant``, ``usertable``, ``userfunc`` and ``usergeometry``.                         |
|                |                         | | ``constant``, ``usertable``, ``userfunc`` assumes that :math:`F_\mathrm{irr}` is the irradiation flux from   |
|                |                         | | from a central star of radius :math:`R_s` and temperature :math:`T_s`, which are set by the second and third |
|                |                         | | parameters respectively. The irradiation flux gets absorbed when propagating radially in a medium with an    |
|                |                         | | absorption opacity :math:`\kappa_\mathrm{irr}`.                                                              |
|                |                         | | For ``constant`` :math:`\kappa_\mathrm{irr}` is constant and set by the fourth parameter.                    |
|                |                         | | For ``usertable`` :math:`\kappa_\mathrm{irr}` is interpolated from a usertable. The third and fourth         |
|                |                         | | parameters set the dimension of the table and the name of the files containing the irradiation opacity.      |
|                |                         | | For ``userfunc``, irradiation opacity is defined from a user defined function. In this case, ``RadSource``   |
|                |                         | | class expects a function to be enrolled with ``Radiation[nFrequencies]::EnrollKappairr(KappairrFunc)``.      |
|                |                         | | For ``usergeometry``, the irradiation source can be entirely prescribed. The user actually prescribes the    |
|                |                         | | source term :math:`-\nabla \cdot \mathbf{F}_\mathrm{irr}` through the enrollment of a function               |
|                |                         | | ``Radiation[nFrequencies]::EnrollIrradiation(IrrFunc)`` in the ``RadSource`` class.                          |
+----------------+-------------------------+----------------------------------------------------------------------------------------------------------------+
| rad_limiter    | string, ...             | | (optional, default is simple) Sets the way the physical condition :math:`f\equiv||F||/E_r\leq 1` is enforced.|
|                |                         | | The choice of ``rad_limiter``is of the most importance in optically thin regions.                            |
|                |                         | | Possible values are: ``simple``, ``flatten`` and ``fpreserving``.                                            |
|                |                         | | ``simple`` ensures :math: `f\leq 1` on the reconstructed values at the cell faces. This can lead to beam     |
|                |                         | | widening in the case of the propagation of a highly collimated radiation field.                              |
|                |                         | | ``flatten`` switch to a flat reconstruction scheme whenever  :math: `f\leq 1` on the faces.                  |
|                |                         | | ``fpreserving`` sets :math: `f\leq 1` by ensuring that the value of :math: `f` is equal on the faces and at  |
|                |                         | | the center of a cell. It is the most accurate implementation.                                                |
+----------------+-------------------------+----------------------------------------------------------------------------------------------------------------+

Using the Radiation module
---------------------

Several examples are provided in the :file:`test/Radiation` directory. Radiation is considered in Idefix as a instance of the `Fluid` class, hence
one can apply the technics used for the gas to the radiation fluid. In order to prepare for future developments, *Idefix* is already able to handle an arbitrarily number of radiation frequency group.
Each frequency group is stored in an instance of `Fluid` and stored in a container (:code:`std::vector radiation`) in the `DataBlock`. The same is true for the mirror `DataBlockHost`: the
radiation primitive variable are all stored in :code:`std::vector RadVc`. However, we emphasize that the current implementation of our source terms only work for one frequency group.

Initialising a single frequency group is done as follow:

.. code-block:: c++


    void Setup::InitFlow(DataBlock &data) {
      // Create a host copy
      DataBlockHost d(data);

      for(int k = 0; k < d.np_tot[KDIR] ; k++) {
          for(int j = 0; j < d.np_tot[JDIR] ; j++) {
              for(int i = 0; i < d.np_tot[IDIR] ; i++) {

                  d.Vc(RHO,k,j,i) = 1.0;            // Set the gas density to 1
                  d.RadVc[0](ER,k,j,i) = 1.0;     // Set first radiation group energy density to 1

                  d.Vc(VX1,k,j,i) = 1;              // Set the gas velocity to 1
                  d.RadVc[0](FR1,k,j,i) = 0.0;     // Set the radiation energy flux to 0

              }
          }
      }

      // Send it all, if needed
      d.SyncToDevice();
    }



All of the radiation fields are automatically outputed in the dump and vtk outputs created by *Idefix*.
