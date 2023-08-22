.. _eosModule:

Equation of state module
=========================

About
---------

By default, *Idefix* can handle either isothermal equation of states (in which case the pressure is never computed, and the code
works with a prescribed sound speed function), or an ideal adiabatic equation of state (assuming a constant adiabatic exponent :math:`\gamma`).

If one wants to compute the dynamics of more complex fluids (e.g. multiphase flows, partial ionisation, etc.), then the ideal adiabatic equation of
state is not sufficient and one needs to code a *custom* equation of state. This is done by implementing the class ``EquationOfState`` with the functions
required by *Idefix* algorithm.

Functions needed
-----------------

*Idefix* requires the definition of 3 functions that are called during the integration loop: a function to compute the first adiabatic exponent :math:`\Gamma_1`, and
the functions needed to convert pressure to internal energy and *vice-versa*. The definition for these function is:

.. code-block:: c++

  class EquationOfState {
    // First adiabatic exponent.
    KOKKOS_INLINE_FUNCTION real GetGamma(real P , real rho ) const {
      real gamma;
      // Compute gamma (needed for sound speed estimations, 5/3 would work in general)
      return gamma;
    }

    // Compute the internal energy from pressure and density
    KOKKOS_INLINE_FUNCTION
    real GetInternalEnergy(real P, real rho) const {
      real eint; // = ...
      return eint;
    }

    // Compute the pressure from internal energy and density
    KOKKOS_INLINE_FUNCTION
    real GetPressure(real Eint, real rho) const {
      real P; // = ...
      return P;
    }
  }


How to use a custom EOS
-----------------------

#. Copy the template file ``eos_template.hpp`` (in src/fluid/eos) in your problem directory and rename it (e.g. ``my_eos.hpp``)
#. Make sure that you have not enabled the ISOTHERMAL approximation in your ``definitions.hpp``
#. Implement your EOS in ``my_eos.hpp``, and in particular the 3 EOS functions required.
#. in cmake, enable ``Idefix_CUSTOM_EOS`` and set ``Idefix_CUSTOM_EOS_FILE`` to ``my_eos.hpp`` (or the filename you have chosen in #1)
#. Compile and run
