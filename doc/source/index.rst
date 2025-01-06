.. Idefix documentation master file, created by
   sphinx-quickstart on Mon Sep 21 10:36:16 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

##############################
Idefix code documentation
##############################
=================
About *Idefix*
=================
*Idefix* is designed to be a performance-portable astrophysical code. This means that it can run both on your laptop's cpu or on the largest GPU HPCs recently
bought by your university. More technically, *Idefix* can run in serial, use OpenMP and/or MPI (message passing interface) for parallelization, and use CUDA or Xeon-Phi for
out of chip acceleration. Of course, all these capabilities are embedded within one single code, so the code relies on relatively abstracted classes and objects available in C++17, which are not necessarily
familiar to astrophysicists. A large effort has been devoted to simplify this level of abstraction so that the code can be modified by researchers and students familiar with C and who are aware of basic object-oriented concepts.

In order to ease the transition from more common codes, *Idefix* has been designed to be close to the PLUTO code. Several of the data structures have identical names, input files are very similar and the numerical integration
algorithm is essentially the same (but with some major modification to its structure). Hence, your favourite PLUTO setup can be converted into an *Idefix* setup within a couple of hours.

================
Requirements
================
*Idefix* is written is standard C++17 and does not rely on any external library in serial (non MPI).

Compiler
  *Idefix* requires a C++17 compatible compiler. It has been tested successfully with GCC (>8), Intel compiler suite (>2018) and
  Clang on both Intel and AMD CPUs. *Idefix* has also been tested on NVIDIA GPUs (Pascal, Volta and Ampere architectures) using the nvcc (>10) compiler, and on AMD GPUs (Radeon Mi50, Mi210, Mi250) using the hipcc compiler.

Kokkos library
  *Idefix* relies internally on the `Kokkos <https://github.com/kokkos/kokkos>`_ library, which is bundled with *Idefix* as a git submodule and compiled on the fly, hence no external installation is required.

MPI library
  When using MPI parallelisation, *Idefix* relies on an external MPI library. *Idefix* has been tested successfully with OpenMPI and IntelMPI libraries. When used on GPU architectures, *Idefix* assumes that
  the MPI library is GPU-Aware. If unsure, check this last point with your system administrator.

Python
  When using *Idefix* with its python interface through the module `Pydefix`, *Idefix* relies on an external python>=3.8 interpreter with the module `pybind11 <https://pybind11.readthedocs.io>`_
  installed.

================
Features
================
* Compressible hydrodynamics in 1D, 2D, 3D
* Compressible magnetohydrodynamics using constrained transport in 1D, 2D, 3D
* Multiple geometry (cartesian, polar, spherical)
* Variable mesh spacing
* Multiple parallelisation strategies (OpenMP, MPI, GPU offloading, etc...)
* Full non-ideal MHD (Ohmic, ambipolar, Hall)
* Viscosity and thermal diffusion
* Super-timestepping for all parabolic terms
* Orbital advection (Fargo-like)
* Self-gravity
* Multi dust species modelled as pressureless fluids
* Multiple planets interraction


===========================
Terms and condition of Use
===========================
*Idefix* is distributed freely under the `CeCILL license <https://en.wikipedia.org/wiki/CeCILL>`_, a free software license adapted to both international and French legal matters, in the spirit of and retaining
compatibility with the GNU General Public License (GPL). We expect *Idefix* to be referenced and acknowledeged by authors in their publications. At the minimum, the authors
should cite the *Idefix* `method paper <https://ui.adsabs.harvard.edu/abs/2023A%26A...677A...9L/abstract>`_.

*Idefix* data structure and algorithm are derived from Andrea Mignone's `PLUTO code <http://plutocode.ph.unito.it/>`_, released under the GPL license.
*Idefix* also relies on the `Kokkos <https://github.com/kokkos/kokkos>`_ performance portability programming ecosystem released under the terms
of Contract DE-NA0003525 with National Technology & Engineering Solutions of Sandia, LLC (NTESS).

==================
Main Contributors
==================

Geoffroy Lesur
  code design and architecture

Soufiane Baghdadi
  code optimisation, 2D Riemann solvers, RKL scheme

Gaylor Wafflard-Fernandez
  planet-disc interaction

Jonah Mauxion
  self-gravity module

Clément Robert
  gitlab integration, linter

Jean Kempf & François Rincon
    anisotropic diffusion

========================
About this documentation
========================

This documentation has automatically been generated on |today| from the following *Idefix* commit:

.. git_commit_detail::
    :branch:
    :commit:
    :sha_length: 10
    :no_github_link:

===================
Acknowledgements
===================

The developement of *Idefix* was supported by the European Research Council (ERC)
under the European Union Horizon 2020 research and innovation programme (Grant agreement No. 815559 (MHDiscs)).
Idefix developement team is partly funded by the `PEPR Origins <https://pepr-origins.fr>`_ through the project "MHD@Exascale".
The Idefix collaboration benefited from funding from the “Programme National de Physique Stellaire” (PNPS),
“Programme National Soleil-Terre” (PNST), “Programme National de Hautes Energies” (PNHE) and
“Programme National de Planétologie” (PNP) of CNRS/INSU co-funded by CEA and CNES.


.. toctree::
   :maxdepth: 3
   :caption: Contents:

   quickstart
   reference
   modules
   programmingguide
   performances
   kokkos
   contributing
   faq
   changelog

   api/library_root




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
