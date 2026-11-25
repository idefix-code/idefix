.. _forcingModule:

Stochastic forcing module
=========================

The forcing module drives turbulence by adding a stochastic external force to the fluid equations.
The force is a superposition of Fourier modes :math:`e^{i\,\mathbf{k}\cdot\mathbf{x}}` whose complex
amplitudes follow independent Ornstein–Uhlenbeck (OU) processes, one per mode and per vector component.
It is typically used to sustain turbulence in decaying-turbulence setups (see e.g. ``test/HD/forcing``).

This approach is identical to the one described in `St-Onge, D. A., M. W. Kunz, J. Squire, et
A. A. Schekochihin, *Journal of Plasma
Physics* 86 (octobre 2020), 905860503 <https://doi.org/10.1017/S0022377820000860>`_.

The forcing term is added to the right-hand side of the momentum equation
(:math:`\partial_t(\rho \mathbf{v}) \mathrel{+}= \rho\,\mathbf{f}`), and to the energy equation
(:math:`\rho\,\mathbf{v}\cdot\mathbf{f}`) when the pressure equation is enabled.

Input parameters
----------------

All parameters are set in the ``[Forcing]`` section of the input file (see :ref:`forcingSection` for
the full list):

``iso3D kmin kmax``
    3D isotropic forcing: all wavevectors with :math:`k_\mathrm{min} \le |\mathbf{k}| \le
    k_\mathrm{max}` are driven. Cartesian geometry with 3 dimensions and 3 components only.

``iso2D normal kmin kmax``
    2D isotropic forcing in the plane perpendicular to ``normal`` (``IDIR``, ``JDIR`` or ``KDIR``),
    with :math:`k_\mathrm{min} \le |\mathbf{k}| \le k_\mathrm{max}`.

``solenoidal``
    When present, the forcing is projected onto its solenoidal (divergence-free) part, computed as
    the curl of the pristine forcing normalised by the mean wavenumber. Not available in polar and
    cylindrical geometries.

``mach M``
    Alternative to ``t_corr``/``epsilon``: the amplitude is set from a target Mach number ``M``,
    using the volume-averaged sound speed :math:`c_s` at :math:`t=0` (``csiso`` in isothermal runs):

    .. math::
        \epsilon = \frac{(M c_s)^3 k_f}{2\pi}, \qquad
        t_\mathrm{corr} = \frac{2\pi}{M c_s k_f}, \qquad
        k_f = \frac{k_\mathrm{min}+k_\mathrm{max}}{2}.

``t_corr`` and ``epsilon``
    OU correlation time and amplitude, used when ``mach`` is not set. ``epsilon`` is internally
    divided by the number of forcing modes, so the total injected amplitude does not depend on how
    many modes are used.

``seed``
    Seed of the random number generator (default: 0).

``stoptime``
    Time after which the forcing is switched off (default: never).

``write``
    If set to 1, the time series of the OU amplitudes and of the normal draws are written to
    ``ou_prank0_seed<seed>.dat`` and ``normal_prank0_seed<seed>.dat`` in the output folder.

How the OU processes are updated
--------------------------------

At each timestep, the complex amplitude :math:`A` of each (mode, component) pair is updated as

.. math::
    A_{n+1} = \frac{\bar{A}}{\sqrt{2}}
        + \left(A_n - \frac{\bar{A}}{\sqrt{2}}\right) e^{-\Delta t / t_\mathrm{corr}}
        + \sqrt{\frac{\epsilon}{\sqrt{2}\,t_\mathrm{corr}}
                 \left(1 - e^{-2\Delta t / t_\mathrm{corr}}\right)}\, \mathcal{N},

where :math:`\bar{A}` is the (zero) mean and :math:`\mathcal{N}` a standard normal draw. The
resulting forcing field at a cell is :math:`\mathbf{f}(\mathbf{x}) =
\sum_l \mathrm{Re}\left[A_l\, e^{i\,\mathbf{k}_l\cdot\mathbf{x}}\right]`.

Random numbers and reproducibility
----------------------------------

Random draws are produced by ``idfx::RandomXorShift64`` (``src/utils/random.hpp``): one
self-managed xorshift64 stream per (mode, component), with standard normal draws obtained through
the Box–Muller transform. Unlike a Kokkos random pool, the sequence is deterministic on every
backend, thread count and MPI rank, and identical on all ranks, and we can restore easily
its state from a restart dump.

The complete forcing state (OU amplitudes ``ouValues``, random state ``ouRngState``, and the OU
parameters ``ouMeans``, ``ouTcorrs``, ``ouEpsilons``) is stored in the restart dumps, so a restart
is bit-exact: the forcing applied after a restart is identical to that of an uninterrupted run.
